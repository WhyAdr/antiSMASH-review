"""Deterministic cohort loading and product/domain matrix construction."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from .clusterblast import ClusterBlastParseError
from .discovery import discover
from .fingerprints import normalize_label, product_counter
from .genbank import GenBankParseError
from .loading import load_review_input
from .models import Record

CohortValueMode = Literal["binary", "count"]


class CohortError(ValueError):
    """Raised when cohort input or member loading is invalid."""


@dataclass(slots=True)
class CohortMember:
    name: str
    input_path: Path
    records: list[Record]
    product_counts: Counter[str]
    domain_counts: Counter[str]
    input_paths: frozenset[Path]


@dataclass(slots=True, frozen=True)
class SkippedCohortMember:
    name: str
    input_path: Path
    error: str


@dataclass(slots=True)
class CohortResult:
    root: Path
    members: list[CohortMember]
    product_columns: list[str]
    domain_columns: list[str]
    product_display_labels: dict[str, str]
    domain_display_labels: dict[str, str]
    product_raw_labels: dict[str, tuple[str, ...]]
    domain_raw_labels: dict[str, tuple[str, ...]]
    product_matrix: list[list[int]]
    domain_matrix: list[list[int]]
    value_mode: CohortValueMode
    input_paths: frozenset[Path] = frozenset()
    skipped: list[SkippedCohortMember] = field(default_factory=list)
    domain_jaccard: list[list[float]] | None = None
    cluster_order: list[str] | None = None
    cluster_newick: str | None = None
    manifest_path: Path | None = None


@dataclass(slots=True, frozen=True)
class _MemberSpec:
    name: str
    input_path: Path


def _domain_label(name: str | None, domain_id: str | None) -> str:
    if name and name.strip():
        return name
    if domain_id and domain_id.strip():
        return domain_id
    return "<unnamed-domain>"


def _labels(records: Iterable[Record], *, domains: bool) -> dict[str, set[str]]:
    labels: dict[str, set[str]] = {}
    for record in records:
        if domains:
            values = (_domain_label(domain.name, domain.domain_id) for domain in record.domains)
        else:
            values = (product for region in record.regions for product in region.products)
        for value in values:
            key = normalize_label(value)
            labels.setdefault(key, set()).add(value)
    return labels


def _display_labels(
    raw_labels: dict[str, set[str]],
) -> tuple[dict[str, str], dict[str, tuple[str, ...]]]:
    display: dict[str, str] = {}
    preserved: dict[str, tuple[str, ...]] = {}
    for key, values in sorted(raw_labels.items()):
        ordered = tuple(sorted(values))
        preserved[key] = ordered
        display[key] = min(
            ordered,
            key=lambda value: (len(value.strip()), normalize_label(value), value),
        ).strip()
    return display, preserved


def _matrix(
    members: list[CohortMember],
    columns: list[str],
    *,
    domains: bool,
    value_mode: CohortValueMode,
) -> list[list[int]]:
    rows: list[list[int]] = []
    for member in members:
        counts = member.domain_counts if domains else member.product_counts
        rows.append(
            [
                counts.get(column, 0) if value_mode == "count" else int(counts.get(column, 0) > 0)
                for column in columns
            ]
        )
    return rows


def _parse_manifest(path: Path) -> list[_MemberSpec]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise CohortError(f"could not read cohort manifest {path}: {exc}") from exc

    specs: list[_MemberSpec] = []
    seen: set[str] = set()
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = line.split("\t")
        if line_number == 1 and [field.strip().casefold() for field in fields] == [
            "sample",
            "path",
        ]:
            continue
        if len(fields) != 2:
            raise CohortError(
                f"invalid cohort manifest line {line_number}: expected sample<TAB>path"
            )
        name, raw_path = (field.strip() for field in fields)
        if not name or not raw_path:
            raise CohortError(f"invalid cohort manifest line {line_number}: empty sample or path")
        key = normalize_label(name)
        if key in seen:
            raise CohortError(f"duplicate cohort sample name: {name}")
        seen.add(key)
        member_path = Path(raw_path)
        if not member_path.is_absolute():
            member_path = path.parent / member_path
        specs.append(_MemberSpec(name=name, input_path=member_path.resolve()))
    if not specs:
        raise CohortError(f"cohort manifest has no members: {path}")
    return specs


def _member_specs(
    root: Path | None,
    manifest_path: Path | None,
) -> tuple[Path, list[_MemberSpec], Path | None]:
    if (root is None) == (manifest_path is None):
        raise CohortError("provide exactly one of cohort root or --manifest")
    if manifest_path is not None:
        manifest = manifest_path.resolve()
        if not manifest.is_file():
            raise CohortError(f"cohort manifest is not a file: {manifest}")
        return manifest.parent, _parse_manifest(manifest), manifest

    assert root is not None
    cohort_root = root.resolve()
    if not cohort_root.is_dir():
        raise CohortError(f"cohort root is not a directory: {cohort_root}")
    children = sorted(
        (path for path in cohort_root.iterdir() if path.is_dir()),
        key=lambda path: (normalize_label(path.name), path.name),
    )
    if not children:
        raise CohortError(f"cohort root has no member directories: {cohort_root}")
    return cohort_root, [_MemberSpec(path.name, path.resolve()) for path in children], None


def _load_member(spec: _MemberSpec, *, lenient: bool) -> tuple[CohortMember, set[Path]]:
    try:
        manifest = discover(spec.input_path, recursive=True)
        loaded = load_review_input(manifest, lenient=lenient)
    except (ClusterBlastParseError, GenBankParseError, OSError, ValueError) as exc:
        raise CohortError(f"member {spec.name!r} at {spec.input_path}: {exc}") from exc

    member = CohortMember(
        name=spec.name,
        input_path=spec.input_path,
        records=loaded.records,
        product_counts=product_counter(loaded.records, normalized=True),
        domain_counts=Counter(
            normalize_label(domain.name or domain.domain_id or "<unnamed-domain>")
            for record in loaded.records
            for domain in record.domains
        ),
        input_paths=frozenset(loaded.input_paths),
    )
    return member, set(loaded.input_paths)


def build_cohort(
    root: Path | None = None,
    *,
    manifest: Path | None = None,
    value_mode: CohortValueMode = "binary",
    lenient: bool = False,
    skip_invalid_members: bool = False,
    cluster_by: Literal["none", "domain-jaccard"] = "none",
) -> CohortResult:
    """Load cohort members and build deterministic product/domain matrices.

    Directory mode sorts immediate member directories. Manifest mode preserves the
    manifest's listed order. Products count regions; domains count adapted aSDomain
    features. Normalized keys use Unicode NFKC, stripping, and case folding while
    retaining raw labels in the result metadata.
    """

    if value_mode not in {"binary", "count"}:
        raise CohortError(f"unknown cohort value mode: {value_mode}")
    if cluster_by not in {"none", "domain-jaccard"}:
        raise CohortError(f"unknown cohort clustering mode: {cluster_by}")
    cohort_root, specs, manifest_path = _member_specs(root, manifest)
    members: list[CohortMember] = []
    skipped: list[SkippedCohortMember] = []
    input_paths: set[Path] = {manifest_path} if manifest_path is not None else set()

    for spec in specs:
        try:
            member, member_input_paths = _load_member(spec, lenient=lenient)
        except CohortError as exc:
            if not skip_invalid_members:
                raise
            skipped.append(SkippedCohortMember(spec.name, spec.input_path, str(exc)))
            continue
        members.append(member)
        input_paths.update(member_input_paths)

    if not members:
        details = f"; skipped {len(skipped)} invalid member(s)" if skipped else ""
        raise CohortError(f"no valid cohort members loaded{details}")

    product_display, product_raw = _display_labels(
        _labels((record for member in members for record in member.records), domains=False)
    )
    domain_display, domain_raw = _display_labels(
        _labels((record for member in members for record in member.records), domains=True)
    )
    product_columns = sorted(product_display)
    domain_columns = sorted(domain_display)

    result = CohortResult(
        root=cohort_root,
        members=members,
        product_columns=product_columns,
        domain_columns=domain_columns,
        product_display_labels=product_display,
        domain_display_labels=domain_display,
        product_raw_labels=product_raw,
        domain_raw_labels=domain_raw,
        product_matrix=_matrix(members, product_columns, domains=False, value_mode=value_mode),
        domain_matrix=_matrix(members, domain_columns, domains=True, value_mode=value_mode),
        value_mode=value_mode,
        input_paths=frozenset(input_paths),
        skipped=skipped,
        manifest_path=manifest_path,
    )
    if cluster_by == "domain-jaccard":
        from .clustering import average_linkage_domain_clustering

        (
            result.domain_jaccard,
            result.cluster_order,
            result.cluster_newick,
        ) = average_linkage_domain_clustering(result.members)
    return result
