from __future__ import annotations

import ast
import hashlib
import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.SeqFeature import CompoundLocation, ExactPosition, SeqFeature
    from Bio.SeqRecord import SeqRecord
except ImportError as exc:  # pragma: no cover - exercised in a clean environment
    raise RuntimeError(
        "GenBank input requires Biopython; install antismash-review with `python -m pip install .`"
    ) from exc

from .locations import overlaps
from .models import (
    AntiSmashProvenance,
    CollectionFeature,
    Diagnostic,
    Domain,
    Gene,
    GeneFunction,
    Location,
    LocationPart,
    Module,
    Motif,
    PfamHit,
    Qualifiers,
    RawFeature,
    Record,
    Severity,
)

_ANTISMASH_DATA_RE = re.compile(
    r"##antiSMASH-Data-START##(?P<body>.*?)##antiSMASH-Data-END##",
    re.IGNORECASE | re.DOTALL,
)
_METADATA_LINE_RE = re.compile(r"^\s*(?P<key>[^:]+?)\s*::\s*(?P<value>.*?)\s*$")
_GENE_FUNCTION_RE = re.compile(
    r"^(?P<category>[\w-]+)(?:\s+\((?P<tool>[^)]+)\))?(?:\s+(?P<description>.*))?$"
)

_ADAPTED_FEATURE_TYPES = frozenset(
    {
        "region",
        "cand_cluster",
        "protocluster",
        "proto_core",
        "CDS",
        "aSDomain",
        "aSModule",
        "CDS_motif",
        "PFAM_domain",
    }
)
_STRUCTURAL_RAW_ONLY_FEATURE_TYPES = frozenset({"source"})


class GenBankParseError(RuntimeError):
    """Raised when a GenBank file cannot be parsed or adapted strictly."""


def _qualifiers(feature: SeqFeature) -> Qualifiers:
    return {
        key: tuple(str(value) for value in values) for key, values in feature.qualifiers.items()
    }


def _values(qualifiers: Mapping[str, Sequence[str]], key: str) -> list[str]:
    return [str(value) for value in qualifiers.get(key, ())]


def _first(
    qualifiers: Mapping[str, Sequence[str]], key: str, default: str | None = None
) -> str | None:
    values = qualifiers.get(key, ())
    return str(values[0]) if values else default


def _first_alias(qualifiers: Mapping[str, Sequence[str]], *keys: str) -> str | None:
    for key in keys:
        value = _first(qualifiers, key)
        if value is not None:
            return value
    return None


def _compact(value: str | None) -> str | None:
    return re.sub(r"\s+", "", value) if value is not None else None


def _compact_values(values: Iterable[str]) -> list[str]:
    return [re.sub(r"\s+", "", value) for value in values]


def _as_int(value: str | None) -> int | None:
    try:
        return int(value) if value is not None else None
    except ValueError:
        return None


def _as_float(value: str | None) -> float | None:
    try:
        return float(value) if value is not None else None
    except ValueError:
        return None


def _as_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    lowered = value.casefold()
    if lowered in {"true", "yes", "1"}:
        return True
    if lowered in {"false", "no", "0"}:
        return False
    return None


def _flag(qualifiers: Mapping[str, Sequence[str]], key: str) -> bool:
    return key in qualifiers


def _location(feature: SeqFeature, record_length: int) -> Location:
    source = feature.location
    if source is None:
        raise ValueError(f"{feature.type} feature has no location")
    parts = tuple(
        LocationPart(
            start=int(part.start),
            end=int(part.end),
            strand=part.strand,
            fuzzy_start=not isinstance(part.start, ExactPosition),
            fuzzy_end=not isinstance(part.end, ExactPosition),
        )
        for part in source.parts
    )
    if not parts:
        raise ValueError(f"{feature.type} feature has an empty location")
    cross_origin = (
        isinstance(source, CompoundLocation)
        and any(part.start == 0 for part in parts)
        and any(part.end == record_length for part in parts)
    )
    return Location(
        start=min(part.start for part in parts),
        end=max(part.end for part in parts),
        strand=source.strand,
        parts=parts,
        cross_origin=cross_origin,
        original=str(source),
    )


def _metadata_value(value: object) -> tuple[str, ...]:
    if isinstance(value, (list, tuple)):
        return tuple(str(item) for item in value)
    if isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            parsed = None
        if isinstance(parsed, (list, tuple)) and all(isinstance(item, str) for item in parsed):
            return tuple(parsed)
    return (str(value),)


def _parse_antismash_metadata_block(comment: str) -> dict[str, tuple[str, ...]]:
    block = _ANTISMASH_DATA_RE.search(comment)
    if block is None:
        return {}
    raw_fields: dict[str, list[str]] = {}
    for line in block.group("body").splitlines():
        match = _METADATA_LINE_RE.match(line)
        if match is None:
            continue
        raw_fields.setdefault(match.group("key").strip(), []).extend(
            _metadata_value(match.group("value").strip())
        )
    return {key: tuple(values) for key, values in raw_fields.items()}


def _raw_antismash_metadata(path: Path) -> list[dict[str, tuple[str, ...]]]:
    """Read raw metadata blocks so repeated keys lost by Biopython survive."""

    text = path.read_text(encoding="utf-8", errors="replace")
    chunks = re.split(r"(?m)^\s*//\s*$", text)
    return [_parse_antismash_metadata_block(chunk) for chunk in chunks if chunk.strip()]


def _antismash_metadata(
    record: SeqRecord,
    raw_fields: Mapping[str, tuple[str, ...]] | None = None,
) -> dict[str, tuple[str, ...]]:
    """Return all antiSMASH structured-comment fields without dropping repeats."""

    structured = record.annotations.get("structured_comment", {})
    if isinstance(structured, Mapping):
        for section_name, section in structured.items():
            if str(section_name).casefold() != "antismash-data":
                continue
            if isinstance(section, Mapping):
                structured_fields = {
                    str(key): _metadata_value(value)
                    for key, value in section.items()
                    if value is not None
                }
                if structured_fields:
                    if raw_fields:
                        return {**structured_fields, **raw_fields}
                    return structured_fields

    if raw_fields:
        return dict(raw_fields)
    return _parse_antismash_metadata_block(str(record.annotations.get("comment", "")))


def _metadata_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.casefold())


def _first_metadata_value(
    fields: Mapping[str, tuple[str, ...]],
    aliases: set[str],
) -> str | None:
    for key, values in fields.items():
        if _metadata_key(key) in aliases and values:
            return values[0]
    return None


def _antismash_provenance(
    record: SeqRecord,
    raw_fields: Mapping[str, tuple[str, ...]] | None = None,
) -> AntiSmashProvenance:
    fields = _antismash_metadata(record, raw_fields)
    database_versions = {
        key: values[0]
        for key, values in fields.items()
        if values and "database" in _metadata_key(key)
    }
    pfam_version = next(
        (values[0] for key, values in fields.items() if values and "pfam" in _metadata_key(key)),
        None,
    )
    detection_rule_set_version = next(
        (
            values[0]
            for key, values in fields.items()
            if values and "detection" in _metadata_key(key) and "rule" in _metadata_key(key)
        ),
        None,
    )
    return AntiSmashProvenance(
        version=_first_metadata_value(fields, {"version", "antismashversion"}),
        run_date=_first_metadata_value(fields, {"rundate"}),
        pfam_version=pfam_version,
        detection_rule_set_version=detection_rule_set_version,
        database_versions=database_versions,
        raw_fields=fields,
    )


def _antismash_version(record: SeqRecord) -> str | None:
    return _antismash_provenance(record).version


def _annotation_text(record: SeqRecord, key: str) -> str | None:
    value = record.annotations.get(key)
    return str(value) if value is not None else None


def _annotation_list(record: SeqRecord, key: str) -> list[str]:
    value = record.annotations.get(key, [])
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value]
    return [str(value)] if value else []


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _numbers(values: Iterable[str]) -> list[int]:
    parsed: list[int] = []
    for value in values:
        number = _as_int(_compact(value))
        if number is not None:
            parsed.append(number)
    return parsed


def _parse_gene_function(raw: str) -> GeneFunction:
    match = _GENE_FUNCTION_RE.match(raw)
    if not match:
        return GeneFunction("unknown", None, None, raw)
    return GeneFunction(
        category=match.group("category"),
        tool=match.group("tool"),
        description=match.group("description"),
        raw=raw,
    )


def _collection(raw: RawFeature) -> CollectionFeature:
    q = raw.qualifiers
    number_keys = {
        "region": "region_number",
        "cand_cluster": "candidate_cluster_number",
        "protocluster": "protocluster_number",
        "proto_core": "protocluster_number",
    }
    reference_keys: dict[str, str | None] = {
        "region": "candidate_cluster_numbers",
        "cand_cluster": "protoclusters",
        "protocluster": "candidate_cluster_numbers",
        "proto_core": None,
    }
    reference_key = reference_keys[raw.feature_type]
    references = _numbers(_values(q, reference_key)) if reference_key else []
    rules = _values(q, "rules") + _values(q, "detection_rules") + _values(q, "detection_rule")
    return CollectionFeature(
        feature_type=raw.feature_type,
        number=_as_int(_compact(_first(q, number_keys[raw.feature_type]))),
        location=raw.location,
        products=_values(q, "product"),
        references=references,
        kind=_first(q, "kind"),
        category=_first(q, "category"),
        rules=rules,
        smiles=_compact_values(_values(q, "SMILES")),
        polymer=_values(q, "polymer"),
        core_location=_compact(_first(q, "core_location")),
        cutoff=_as_int(_first(q, "cutoff")),
        neighbourhood=_as_int(_first(q, "neighbourhood")),
        creating_tool=_first_alias(q, "aSTool", "aStool"),
        contig_edge=_as_bool(_first(q, "contig_edge")),
        qualifiers=q,
    )


def _gene(raw: RawFeature) -> Gene:
    q = raw.qualifiers
    return Gene(
        location=raw.location,
        locus_tag=_first(q, "locus_tag"),
        gene=_first(q, "gene"),
        product=_first(q, "product"),
        protein_id=_first(q, "protein_id"),
        translation=_compact(_first(q, "translation")),
        gene_kind=_first(q, "gene_kind", "unclassified") or "unclassified",
        gene_functions=[_parse_gene_function(value) for value in _values(q, "gene_functions")],
        ec_numbers=_values(q, "EC_number"),
        db_xrefs=_values(q, "db_xref"),
        notes=_values(q, "note"),
        inference=_values(q, "inference"),
        region_numbers=[],
        candidate_cluster_numbers=[],
        protocluster_numbers=[],
        proto_core_numbers=[],
        qualifiers=q,
    )


def _domain(raw: RawFeature) -> Domain:
    q = raw.qualifiers
    subtypes = _values(q, "domain_subtypes") or _values(q, "domain_subtype")
    return Domain(
        domain_id=_compact(_first(q, "domain_id")),
        name=_first(q, "aSDomain"),
        subtypes=subtypes,
        tool=_first_alias(q, "aSTool", "aStool"),
        locus_tag=_first(q, "locus_tag"),
        score=_as_float(_first(q, "score")),
        evalue=_as_float(_first(q, "evalue")),
        protein_start=_as_int(_first(q, "protein_start")),
        protein_end=_as_int(_first(q, "protein_end")),
        specificity=_values(q, "specificity"),
        location=raw.location,
        qualifiers=q,
    )


def _module(raw: RawFeature) -> Module:
    q = raw.qualifiers
    return Module(
        location=raw.location,
        domain_ids=_compact_values(_values(q, "domains")),
        locus_tags=_compact_values(_values(q, "locus_tags")),
        module_type=_first(q, "type"),
        complete=_flag(q, "complete") and not _flag(q, "incomplete"),
        starter=_flag(q, "starter_module"),
        final=_flag(q, "final_module"),
        iterative=_flag(q, "iterative"),
        monomer_pairings=_values(q, "monomer_pairings"),
        multi_cds=len(_values(q, "locus_tags")) > 1,
        qualifiers=q,
    )


def _motif(raw: RawFeature) -> Motif:
    q = raw.qualifiers
    return Motif(
        location=raw.location,
        label=_first(q, "motif") or _first(q, "label"),
        locus_tag=_first(q, "locus_tag"),
        tool=_first_alias(q, "aSTool", "aStool"),
        score=_as_float(_first(q, "score")),
        evalue=_as_float(_first(q, "evalue") or _first(q, "evalue_label")),
        prepeptide=_first(q, "prepeptide"),
        core_sequence=_compact(_first(q, "core_sequence")),
        qualifiers=q,
    )


def _pfam_hit(raw: RawFeature) -> PfamHit:
    q = raw.qualifiers
    return PfamHit(
        location=raw.location,
        accession=_first(q, "db_xref"),
        description=_first(q, "description"),
        locus_tag=_first(q, "locus_tag"),
        tool=_first_alias(q, "aSTool", "aStool"),
        score=_as_float(_first(q, "score")),
        evalue=_as_float(_first(q, "evalue")),
        protein_start=_as_int(_first(q, "protein_start")),
        protein_end=_as_int(_first(q, "protein_end")),
        qualifiers=q,
    )


def _adapt_feature(record: Record, raw: RawFeature) -> None:
    if raw.feature_type in {"region", "cand_cluster", "protocluster", "proto_core"}:
        target = {
            "region": record.regions,
            "cand_cluster": record.candidate_clusters,
            "protocluster": record.protoclusters,
            "proto_core": record.proto_cores,
        }[raw.feature_type]
        target.append(_collection(raw))
    elif raw.feature_type == "CDS":
        record.genes.append(_gene(raw))
    elif raw.feature_type == "aSDomain":
        record.domains.append(_domain(raw))
    elif raw.feature_type == "aSModule":
        record.modules.append(_module(raw))
    elif raw.feature_type == "CDS_motif":
        record.motifs.append(_motif(raw))
    elif raw.feature_type == "PFAM_domain":
        record.pfam_hits.append(_pfam_hit(raw))


def _diagnose_unrecognized_features(record: Record) -> None:
    unknown = sorted(
        {
            raw.feature_type
            for raw in record.raw_features
            if raw.feature_type not in _ADAPTED_FEATURE_TYPES
            and raw.feature_type not in _STRUCTURAL_RAW_ONLY_FEATURE_TYPES
        }
    )
    if unknown:
        record.diagnostics.append(
            Diagnostic(
                code="unrecognized_feature_type",
                severity=Severity.NOTICE,
                message=(
                    "Feature types retained only as raw evidence and not adapted: "
                    + ", ".join(unknown)
                ),
                source=str(record.source_path),
                record_id=record.record_id,
            )
        )


def _resolve_modules(record: Record) -> None:
    by_id: dict[str, list[Domain]] = {}
    for domain in record.domains:
        if domain.domain_id:
            by_id.setdefault(domain.domain_id, []).append(domain)

    for module in record.modules:
        module.missing_domain_ids = [
            domain_id for domain_id in module.domain_ids if domain_id not in by_id
        ]
        for domain_id in module.domain_ids:
            matches = by_id.get(domain_id, [])
            if not matches:
                record.diagnostics.append(
                    Diagnostic(
                        code="module_domain_missing",
                        severity=Severity.WARNING,
                        message=f"Module references missing domain {domain_id}",
                        source=str(record.source_path),
                        record_id=record.record_id,
                    )
                )
            elif len(matches) > 1:
                record.diagnostics.append(
                    Diagnostic(
                        code="domain_id_duplicated",
                        severity=Severity.WARNING,
                        message=f"Domain ID occurs {len(matches)} times: {domain_id}",
                        source=str(record.source_path),
                        record_id=record.record_id,
                    )
                )


def _member_numbers(gene: Gene, collections: list[CollectionFeature]) -> list[int]:
    return [
        collection.number
        for collection in collections
        if collection.number is not None and overlaps(collection.location, gene.location)
    ]


def _assign_gene_memberships(record: Record) -> None:
    for gene in record.genes:
        gene.region_numbers = _member_numbers(gene, record.regions)
        gene.candidate_cluster_numbers = _member_numbers(gene, record.candidate_clusters)
        gene.protocluster_numbers = _member_numbers(gene, record.protoclusters)
        gene.proto_core_numbers = _member_numbers(gene, record.proto_cores)


def _parse_record(
    source: SeqRecord,
    path: Path,
    *,
    source_sha256: str,
    lenient: bool,
    raw_metadata: Mapping[str, tuple[str, ...]] | None = None,
) -> Record:
    provenance = _antismash_provenance(source, raw_metadata)
    record = Record(
        record_id=source.id,
        name=source.name,
        description=source.description,
        length=len(source.seq),
        molecule_type=_annotation_text(source, "molecule_type"),
        topology=_annotation_text(source, "topology"),
        source_path=path.resolve(),
        source_sha256=source_sha256,
        antismash_version=provenance.version,
        organism=_annotation_text(source, "organism"),
        taxonomy=_annotation_list(source, "taxonomy"),
        antismash_provenance=provenance,
    )

    for index, feature in enumerate(source.features):
        qualifiers = _qualifiers(feature)
        raw = RawFeature(
            feature_type=feature.type,
            location=_location(feature, len(source.seq)),
            qualifiers=qualifiers,
            feature_index=index,
        )
        record.raw_features.append(raw)
        try:
            _adapt_feature(record, raw)
        except (KeyError, TypeError, ValueError) as exc:
            diagnostic = Diagnostic(
                code="feature_adapter_failed",
                severity=Severity.WARNING if lenient else Severity.ERROR,
                message=f"Could not adapt {feature.type}: {exc}",
                source=str(path),
                record_id=source.id,
                feature_index=index,
            )
            record.diagnostics.append(diagnostic)
            if not lenient:
                raise GenBankParseError(diagnostic.message) from exc

    _resolve_modules(record)
    _assign_gene_memberships(record)
    _diagnose_unrecognized_features(record)
    return record


def parse_genbank(path: Path, *, lenient: bool = False) -> list[Record]:
    """Parse every GenBank record independently using Biopython."""
    path = Path(path)
    records: list[Record] = []
    try:
        source_sha256 = _sha256(path)
        raw_metadata = _raw_antismash_metadata(path)
        parsed_records = SeqIO.parse(path, "genbank")
        for index, source_record in enumerate(parsed_records):
            records.append(
                _parse_record(
                    source_record,
                    path,
                    source_sha256=source_sha256,
                    lenient=lenient,
                    raw_metadata=raw_metadata[index] if index < len(raw_metadata) else None,
                )
            )
    except (OSError, ValueError, TypeError) as exc:
        raise GenBankParseError(f"Could not parse {path}: {exc}") from exc

    if not records:
        raise GenBankParseError(f"No GenBank records found in {path}")
    return records
