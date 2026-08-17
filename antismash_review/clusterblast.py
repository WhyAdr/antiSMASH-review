from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Sequence
from pathlib import Path
from typing import cast

from antismash_review.models import (
    ClusterBlastHit,
    ClusterBlastPairing,
    ClusterBlastResult,
    ClusterBlastSearchType,
    Diagnostic,
    Record,
    Severity,
)

_REGION_NUM_RE = re.compile(r"_c(?P<region>\d+)\.txt$", re.IGNORECASE)
_HEADER_RE = re.compile(r"^ClusterBlast scores for\s+(?P<record_id>\S+)", re.IGNORECASE)
_SIG_HIT_RE = re.compile(r"^(?P<rank>\d+)\.\s+(?P<accession>\S+)(?:\t(?P<description>.*))?$")
_DETAILS_HEADER_RE = re.compile(r"^(?P<rank>\d+)\.\s+(?P<accession>\S+)$")

# Upstream serializer version → antiSMASH generation (documentation/validation).
_CLUSTERBLAST_RESULT_SCHEMA_GENERATIONS: dict[int, str] = {
    1: "antiSMASH 5.x/6.x",
    2: "antiSMASH 7.0.x",
    3: "antiSMASH 7.1.x",
    5: "antiSMASH 8.x",
}


class ClusterBlastParseError(RuntimeError):
    """A ClusterBlast sidecar was recognized but could not be parsed safely."""


def _required_string(value: object, field: str, path: Path) -> str:
    if not isinstance(value, str):
        raise ClusterBlastParseError(f"ClusterBlast {field} is missing or not a string in {path}")
    return value


def _required_nonempty_string(value: object, field: str, path: Path) -> str:
    result = _required_string(value, field, path)
    if not result:
        raise ClusterBlastParseError(f"ClusterBlast {field} is empty in {path}")
    return result


def _optional_string(value: object, field: str, path: Path) -> str | None:
    if value is None:
        return None
    return _required_string(value, field, path)


def _required_integer(value: object, field: str, path: Path) -> int:
    if type(value) is not int:
        raise ClusterBlastParseError(f"ClusterBlast {field} is not an integer in {path}: {value!r}")
    return value


def _optional_integer(value: object, field: str, path: Path) -> int | None:
    if value is None:
        return None
    return _required_integer(value, field, path)


def _optional_nonnegative_integer(value: object, field: str, path: Path) -> int | None:
    result = _optional_integer(value, field, path)
    if result is not None and result < 0:
        raise ClusterBlastParseError(
            f"ClusterBlast {field} must not be negative in {path}: {result}"
        )
    return result


def _required_float(value: object, field: str, path: Path) -> float:
    if type(value) not in {int, float}:
        raise ClusterBlastParseError(f"ClusterBlast {field} is not numeric in {path}: {value!r}")
    try:
        result = float(cast(int | float, value))
    except OverflowError as exc:
        raise ClusterBlastParseError(
            f"ClusterBlast {field} is not finite in {path}: {value!r}"
        ) from exc
    if not math.isfinite(result):
        raise ClusterBlastParseError(f"ClusterBlast {field} is not finite in {path}: {value!r}")
    return result


def _optional_float(value: object, field: str, path: Path) -> float | None:
    if value is None:
        return None
    return _required_float(value, field, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_clusterblast_text(
    path: Path,
    *,
    search_type: ClusterBlastSearchType,
) -> ClusterBlastResult:
    """Parse one antiSMASH ClusterBlast text result."""
    path = Path(path)
    match_region = _REGION_NUM_RE.search(path.name)
    if not match_region:
        raise ClusterBlastParseError(f"Could not extract region number from filename {path.name}")
    region_number = int(match_region.group("region"))

    try:
        content = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ClusterBlastParseError(f"Could not read {path}: {exc}") from exc

    lines = content.splitlines()
    if not lines:
        raise ClusterBlastParseError(f"Empty ClusterBlast file: {path}")

    header_match = _HEADER_RE.match(lines[0].strip())
    if not header_match:
        raise ClusterBlastParseError(
            f"Missing or invalid ClusterBlast header in {path}: {lines[0]!r}"
        )
    record_id = header_match.group("record_id")

    # Locate "Significant hits:" and "Details:"
    sig_index = -1
    details_index = -1
    for idx, raw_line in enumerate(lines):
        line = raw_line.strip()
        if line == "Significant hits:":
            sig_index = idx
        elif line == "Details:":
            details_index = idx

    if sig_index == -1 or details_index == -1 or details_index < sig_index:
        raise ClusterBlastParseError(
            f"Malformed structure in {path}: Significant hits or Details section missing"
        )

    # Parse Significant hits
    sig_hits: list[tuple[int, str, str]] = []
    expected_rank = 1
    for raw_line in lines[sig_index + 1 : details_index]:
        line = raw_line.strip()
        if not line:
            continue
        sig_match = _SIG_HIT_RE.match(line)
        if not sig_match:
            raise ClusterBlastParseError(f"Malformed significant hit line in {path}: {line!r}")
        rank = int(sig_match.group("rank"))
        if rank != expected_rank:
            raise ClusterBlastParseError(
                f"Non-consecutive rank in {path}: expected {expected_rank}, got {rank}"
            )
        accession = sig_match.group("accession")
        description = sig_match.group("description") or ""
        sig_hits.append((rank, accession, description))
        expected_rank += 1

    # Parse Details section
    details_lines = lines[details_index + 1 :]
    # Split details into blocks separated by ">>"
    details_blocks: list[list[str]] = []
    current_block: list[str] = []
    for raw_line in details_lines:
        line = raw_line.rstrip()
        if line.strip() == ">>":
            if current_block:
                details_blocks.append(current_block)
            current_block = []
        else:
            current_block.append(line)
    if current_block:
        details_blocks.append(current_block)

    # Filter out empty blocks
    details_blocks = [b for b in details_blocks if any(line.strip() for line in b)]

    if len(details_blocks) != len(sig_hits):
        raise ClusterBlastParseError(
            f"Mismatched hit count in {path}: {len(sig_hits)} significant hits, "
            f"{len(details_blocks)} detail blocks"
        )

    rankings: list[ClusterBlastHit] = []
    for (sig_rank, sig_accession, sig_description), block in zip(
        sig_hits, details_blocks, strict=True
    ):
        # First non-empty line of block must be "<rank>. <accession>"
        non_empty = [line.strip() for line in block if line.strip()]
        if not non_empty:
            raise ClusterBlastParseError(f"Empty detail block in {path}")
        det_header = _DETAILS_HEADER_RE.match(non_empty[0])
        if not det_header:
            raise ClusterBlastParseError(f"Malformed detail header in {path}: {non_empty[0]!r}")
        det_rank = int(det_header.group("rank"))
        det_accession = det_header.group("accession")
        if det_rank != sig_rank or det_accession != sig_accession:
            raise ClusterBlastParseError(
                f"Detail header {det_rank}. {det_accession} does not match "
                f"significant hit {sig_rank}. {sig_accession} in {path}"
            )

        cluster_type: str | None = None
        num_hits: int | None = None
        blast_score: float | None = None
        pairings: list[ClusterBlastPairing] = []
        in_blast_table = False

        for raw_line in block:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("Type:"):
                cluster_type = line.split(":", 1)[1].strip() or None
            elif line.startswith("Number of proteins with BLAST hits to this cluster:"):
                try:
                    num_hits = int(line.split(":", 1)[1].strip())
                except ValueError as exc:
                    raise ClusterBlastParseError(f"Invalid hit count in {path}: {line!r}") from exc
            elif line.startswith("Cumulative BLAST score:"):
                try:
                    blast_score = float(line.split(":", 1)[1].strip())
                except ValueError as exc:
                    raise ClusterBlastParseError(
                        f"Invalid blast score in {path}: {line!r}"
                    ) from exc
            elif line.startswith("Table of Blast hits"):
                in_blast_table = True
            elif in_blast_table:
                cols = raw_line.split("\t")
                if len(cols) == 6:
                    try:
                        pairings.append(
                            ClusterBlastPairing(
                                query_gene=cols[0].strip(),
                                subject_gene=cols[1].strip(),
                                percent_identity=float(cols[2]),
                                blast_score=float(cols[3]),
                                percent_coverage=float(cols[4]),
                                evalue=float(cols[5]),
                            )
                        )
                    except ValueError as exc:
                        raise ClusterBlastParseError(
                            f"Invalid numeric value in BLAST hit row in {path}: {raw_line!r}"
                        ) from exc
                elif cols and cols[0].strip():
                    raise ClusterBlastParseError(
                        f"Expected 6 tab-separated BLAST hit columns in {path}, "
                        f"got {len(cols)}: {raw_line!r}"
                    )

        rankings.append(
            ClusterBlastHit(
                rank=sig_rank,
                accession=sig_accession,
                description=sig_description,
                cluster_type=cluster_type,
                num_hits=num_hits,
                core_gene_hits=None,
                blast_score=blast_score,
                synteny_score=None,
                core_bonus=None,
                similarity=None,
                pairings=pairings,
            )
        )

    return ClusterBlastResult(
        record_id=record_id,
        region_number=region_number,
        search_type=search_type,
        total_hits=None,
        rankings=rankings,
        source_path=path.resolve(),
        source_sha256=_sha256(path),
        source_format="text",
    )


def _parse_clusterblast_json_unchecked(path: Path) -> list[ClusterBlastResult]:
    """Parse only ClusterBlast modules from a native antiSMASH JSON sidecar."""
    path = Path(path)
    try:
        content = path.read_text(encoding="utf-8")
        data = json.loads(content)
    except (OSError, json.JSONDecodeError) as exc:
        raise ClusterBlastParseError(f"Could not parse JSON {path}: {exc}") from exc

    if not isinstance(data, dict) or "records" not in data:
        raise ClusterBlastParseError(f"Invalid antiSMASH JSON root in {path}")

    raw_records = data["records"]
    if not isinstance(raw_records, list):
        raise ClusterBlastParseError(f"antiSMASH JSON records is not a list in {path}")

    source_sha256 = _sha256(path)
    results: list[ClusterBlastResult] = []

    type_mapping: tuple[tuple[str, ClusterBlastSearchType], ...] = (
        ("general", "clusterblast"),
        ("knowncluster", "knownclusterblast"),
        ("subcluster", "subclusterblast"),
    )

    for record_idx, rec in enumerate(raw_records):
        if not isinstance(rec, dict):
            raise ClusterBlastParseError(
                f"antiSMASH JSON record {record_idx} is not an object in {path}"
            )
        rec_id = _required_nonempty_string(rec.get("id"), "containing record id", path)
        modules = rec.get("modules", {})
        if not isinstance(modules, dict):
            raise ClusterBlastParseError(
                f"antiSMASH JSON modules for record {record_idx} is not an object in {path}"
            )
        if "antismash.modules.clusterblast" not in modules:
            continue

        cb_mod = modules["antismash.modules.clusterblast"]
        if not isinstance(cb_mod, dict):
            raise ClusterBlastParseError(
                f"ClusterBlast module in record {record_idx} is not a dict in {path}"
            )

        module_schema_version = _required_integer(
            cb_mod.get("schema_version"), "module schema_version", path
        )
        if module_schema_version not in {1, 2}:
            raise ClusterBlastParseError(
                f"Unsupported ClusterBlast module schema version {module_schema_version!r} "
                f"(expected 1 or 2) in {path}"
            )

        cb_record_id = _required_nonempty_string(cb_mod.get("record_id"), "module record_id", path)
        if cb_record_id != rec_id:
            raise ClusterBlastParseError(
                f"ClusterBlast module record_id {cb_record_id!r} does not match "
                f"record id {rec_id!r} in {path}"
            )

        for sec_key, search_type in type_mapping:
            if sec_key not in cb_mod:
                continue
            section = cb_mod[sec_key]
            if not isinstance(section, dict):
                raise ClusterBlastParseError(
                    f"ClusterBlast {sec_key} section is not an object in {path}"
                )

            result_schema_version = _required_integer(
                section.get("schema_version"), f"{sec_key} result schema_version", path
            )
            if result_schema_version not in _CLUSTERBLAST_RESULT_SCHEMA_GENERATIONS:
                raise ClusterBlastParseError(
                    f"Unsupported ClusterBlast {sec_key} result schema version "
                    f"{result_schema_version!r} "
                    f"(expected {sorted(_CLUSTERBLAST_RESULT_SCHEMA_GENERATIONS)}) in {path}"
                )

            data_version = _optional_string(
                section.get("data_version"), f"{sec_key} data_version", path
            )

            # Validate optional section-level record_id when present.
            if "record_id" in section and section["record_id"] is not None:
                sec_rec_id = _required_nonempty_string(
                    section["record_id"], f"{sec_key} record_id", path
                )
                if sec_rec_id != cb_record_id:
                    raise ClusterBlastParseError(
                        f"ClusterBlast {sec_key} record_id {sec_rec_id!r} does not match "
                        f"module record_id {cb_record_id!r} in {path}"
                    )

            # Upstream section keys omit the ``blast`` suffix used by the
            # serialized GeneralResults.search_type value (for example,
            # ``general`` contains ``search_type=clusterblast``).
            if "search_type" in section and section["search_type"] is not None:
                sec_search_type = _required_nonempty_string(
                    section["search_type"], f"{sec_key} search_type", path
                )
                if sec_search_type != search_type:
                    raise ClusterBlastParseError(
                        f"ClusterBlast {sec_key} search_type {sec_search_type!r} does not "
                        f"match expected value {search_type!r} in {path}"
                    )

            raw_results = section.get("results")
            if not isinstance(raw_results, list):
                raise ClusterBlastParseError(
                    f"ClusterBlast {sec_key} results is not a list in {path}"
                )

            for region_res in raw_results:
                if not isinstance(region_res, dict):
                    raise ClusterBlastParseError(
                        f"ClusterBlast {sec_key} region result is not an object in {path}"
                    )
                if "region_number" not in region_res:
                    raise ClusterBlastParseError(
                        f"ClusterBlast region result missing region_number in {path}"
                    )
                region_number = _required_integer(
                    region_res["region_number"], "region_number", path
                )
                if region_number < 1:
                    raise ClusterBlastParseError(
                        f"ClusterBlast region_number must be positive in {path}: {region_number}"
                    )
                total_hits = _optional_nonnegative_integer(
                    region_res.get("total_hits"), "total_hits", path
                )
                raw_rankings = region_res.get("ranking", [])
                if not isinstance(raw_rankings, list):
                    raise ClusterBlastParseError(f"ClusterBlast ranking is not a list in {path}")

                rankings: list[ClusterBlastHit] = []
                for rank_idx, rank_entry in enumerate(raw_rankings, start=1):
                    if not isinstance(rank_entry, (list, tuple)) or len(rank_entry) != 2:
                        raise ClusterBlastParseError(
                            f"Invalid ranking entry format in {path}: {rank_entry!r}"
                        )
                    hit_info, hit_details = rank_entry
                    if not isinstance(hit_info, dict) or not isinstance(hit_details, dict):
                        raise ClusterBlastParseError(
                            f"Invalid ranking elements in {path}: {rank_entry!r}"
                        )

                    pairings: list[ClusterBlastPairing] = []
                    raw_pairings = hit_details.get("pairings", [])
                    if not isinstance(raw_pairings, list):
                        raise ClusterBlastParseError(
                            f"ClusterBlast pairings is not a list in {path}"
                        )
                    for pairing_entry in raw_pairings:
                        if not isinstance(pairing_entry, (list, tuple)) or len(pairing_entry) != 3:
                            raise ClusterBlastParseError(
                                f"Invalid pairing entry in {path}: {pairing_entry!r}"
                            )
                        query_str, subject_idx, pairing_dict = pairing_entry
                        if not isinstance(query_str, str):
                            raise ClusterBlastParseError(
                                f"Invalid query string type in {path}: {query_str!r}"
                            )
                        parts = query_str.split("|", 5)
                        if len(parts) >= 5 and parts[4].strip():
                            query_gene = parts[4].strip()
                        elif "|" not in query_str and query_str.strip():
                            query_gene = query_str.strip()
                        else:
                            raise ClusterBlastParseError(
                                f"Malformed query string in {path}: {query_str!r}"
                            )
                        if not isinstance(pairing_dict, dict):
                            raise ClusterBlastParseError(
                                f"Invalid pairing dictionary in {path}: {pairing_dict!r}"
                            )
                        pairings.append(
                            ClusterBlastPairing(
                                query_gene=query_gene,
                                subject_gene=_required_nonempty_string(
                                    pairing_dict.get("name"), "pairing name", path
                                ),
                                percent_identity=_required_float(
                                    pairing_dict.get("perc_ident"), "pairing perc_ident", path
                                ),
                                blast_score=_required_float(
                                    pairing_dict.get("blastscore"), "pairing blastscore", path
                                ),
                                percent_coverage=_required_float(
                                    pairing_dict.get("perc_coverage"),
                                    "pairing perc_coverage",
                                    path,
                                ),
                                evalue=_required_float(
                                    pairing_dict.get("evalue"), "pairing evalue", path
                                ),
                                subject_protein_id=_optional_string(
                                    pairing_dict.get("locus_tag"),
                                    "pairing locus_tag",
                                    path,
                                ),
                                subject_index=_optional_integer(
                                    subject_idx, "pairing subject_index", path
                                ),
                            )
                        )

                    rankings.append(
                        ClusterBlastHit(
                            rank=rank_idx,
                            accession=_required_nonempty_string(
                                hit_info.get("accession"), "hit accession", path
                            ),
                            description=_required_string(
                                hit_info.get("description"), "hit description", path
                            ),
                            cluster_type=_optional_string(
                                hit_info.get("cluster_type"), "hit cluster_type", path
                            ),
                            num_hits=_optional_nonnegative_integer(
                                hit_details.get("hits"), "hit count", path
                            ),
                            core_gene_hits=_optional_nonnegative_integer(
                                hit_details.get("core_gene_hits"),
                                "core gene hit count",
                                path,
                            ),
                            blast_score=_optional_float(
                                hit_details.get("blast_score"), "hit blast_score", path
                            ),
                            synteny_score=_optional_integer(
                                hit_details.get("synteny_score"), "synteny_score", path
                            ),
                            core_bonus=_optional_integer(
                                hit_details.get("core_bonus"), "core_bonus", path
                            ),
                            similarity=_optional_nonnegative_integer(
                                hit_details.get("similarity"), "similarity", path
                            ),
                            pairings=pairings,
                        )
                    )

                results.append(
                    ClusterBlastResult(
                        record_id=cb_record_id,
                        region_number=region_number,
                        search_type=search_type,
                        total_hits=total_hits,
                        rankings=rankings,
                        source_path=path.resolve(),
                        source_sha256=source_sha256,
                        source_format="json",
                        module_schema_version=module_schema_version,
                        result_schema_version=result_schema_version,
                        data_version=data_version,
                    )
                )

    return results


def parse_clusterblast_json(path: Path) -> list[ClusterBlastResult]:
    """Parse supported ClusterBlast JSON or raise one public enrichment error."""
    path = Path(path)
    try:
        return _parse_clusterblast_json_unchecked(path)
    except ClusterBlastParseError:
        raise
    except (KeyError, TypeError, ValueError, AttributeError, OverflowError) as exc:
        raise ClusterBlastParseError(
            f"Malformed ClusterBlast JSON structure in {path}: {exc}"
        ) from exc


def _result_key(
    result: ClusterBlastResult,
) -> tuple[str, int, ClusterBlastSearchType]:
    return (result.record_id, result.region_number, result.search_type)


def merge_clusterblast_results(
    text_results: Sequence[ClusterBlastResult],
    json_results: Sequence[ClusterBlastResult],
    *,
    lenient: bool = False,
    diagnostics: list[Diagnostic] | None = None,
) -> list[ClusterBlastResult]:
    selected: dict[
        tuple[str, int, ClusterBlastSearchType],
        ClusterBlastResult,
    ] = {}
    json_seen: set[tuple[str, int, ClusterBlastSearchType]] = set()
    text_seen: set[tuple[str, int, ClusterBlastSearchType]] = set()

    for result in json_results:
        key = _result_key(result)
        if key in json_seen:
            msg = f"duplicate JSON ClusterBlast result: {key}"
            if not lenient:
                raise ClusterBlastParseError(msg)
            if diagnostics is not None:
                diagnostics.append(
                    Diagnostic(
                        code="clusterblast_duplicate_result",
                        severity=Severity.WARNING,
                        message=msg,
                        source=str(result.source_path),
                        record_id=result.record_id,
                    )
                )
            continue
        json_seen.add(key)
        selected[key] = result

    for result in text_results:
        key = _result_key(result)
        if key in text_seen:
            msg = f"duplicate text ClusterBlast result: {key}"
            if not lenient:
                raise ClusterBlastParseError(msg)
            if diagnostics is not None:
                diagnostics.append(
                    Diagnostic(
                        code="clusterblast_duplicate_result",
                        severity=Severity.WARNING,
                        message=msg,
                        source=str(result.source_path),
                        record_id=result.record_id,
                    )
                )
            continue
        text_seen.add(key)
        selected[key] = result

    return [selected[key] for key in sorted(selected)]


def attach_clusterblast_results(
    records: list[Record],
    results: list[ClusterBlastResult],
    *,
    lenient: bool = False,
    diagnostics: list[Diagnostic] | None = None,
) -> None:
    valid_attachments: list[tuple[Record, ClusterBlastResult]] = []
    for result in results:
        candidates = [
            record
            for record in records
            if result.record_id in {record.record_id, record.name}
            and any(region.number == result.region_number for region in record.regions)
        ]
        if len(candidates) != 1:
            msg = f"expected one GenBank target for {_result_key(result)}, found {len(candidates)}"
            if not lenient:
                raise ClusterBlastParseError(msg)
            diag = Diagnostic(
                code="clusterblast_attach_failed",
                severity=Severity.WARNING,
                message=msg,
                source=str(result.source_path),
                record_id=result.record_id,
            )
            # Route to the owning record if uniquely identifiable.
            id_matches = [r for r in records if result.record_id in {r.record_id, r.name}]
            if len(id_matches) == 1:
                id_matches[0].diagnostics.append(diag)
            elif diagnostics is not None:
                diagnostics.append(diag)
        else:
            valid_attachments.append((candidates[0], result))

    for target_record, result in valid_attachments:
        target_record.clusterblast_results.append(result)
