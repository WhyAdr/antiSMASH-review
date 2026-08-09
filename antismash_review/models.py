from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Literal

ClusterBlastSearchType = Literal[
    "clusterblast",
    "knownclusterblast",
    "subclusterblast",
]
ClusterBlastSourceFormat = Literal["text", "json"]


@dataclass(slots=True, frozen=True)
class ClusterBlastPairing:
    query_gene: str
    subject_gene: str
    percent_identity: float
    blast_score: float
    percent_coverage: float
    evalue: float
    subject_protein_id: str | None = None
    subject_index: int | None = None


@dataclass(slots=True)
class ClusterBlastHit:
    rank: int
    accession: str
    description: str
    cluster_type: str | None
    num_hits: int | None
    core_gene_hits: int | None
    blast_score: float | None
    synteny_score: int | None
    core_bonus: int | None
    similarity: int | None
    pairings: list[ClusterBlastPairing] = field(default_factory=list)


@dataclass(slots=True)
class ClusterBlastResult:
    record_id: str
    region_number: int
    search_type: ClusterBlastSearchType
    total_hits: int | None
    rankings: list[ClusterBlastHit]
    source_path: Path
    source_sha256: str
    source_format: ClusterBlastSourceFormat
    module_schema_version: int | None = None
    result_schema_version: int | None = None


class Severity(str, Enum):
    ERROR = "error"
    WARNING = "warning"
    NOTICE = "notice"


@dataclass(slots=True, frozen=True)
class Diagnostic:
    code: str
    severity: Severity
    message: str
    source: str
    record_id: str | None = None
    feature_index: int | None = None


@dataclass(slots=True, frozen=True)
class LocationPart:
    start: int
    end: int
    strand: int | None
    fuzzy_start: bool = False
    fuzzy_end: bool = False


@dataclass(slots=True, frozen=True)
class Location:
    start: int
    end: int
    strand: int | None
    parts: tuple[LocationPart, ...]
    cross_origin: bool
    original: str

    @property
    def display(self) -> str:
        rendered = ",".join(f"{part.start + 1}..{part.end}" for part in self.parts)
        rendered = f"join({rendered})" if len(self.parts) > 1 else rendered
        return f"complement({rendered})" if self.strand == -1 else rendered

    @property
    def partial(self) -> bool:
        return any(part.fuzzy_start or part.fuzzy_end for part in self.parts)


Qualifiers = dict[str, tuple[str, ...]]


@dataclass(slots=True, frozen=True)
class RawFeature:
    feature_type: str
    location: Location
    qualifiers: Qualifiers
    feature_index: int


@dataclass(slots=True, frozen=True)
class GeneFunction:
    category: str
    tool: str | None
    description: str | None
    raw: str


@dataclass(slots=True)
class Gene:
    location: Location
    locus_tag: str | None
    gene: str | None
    product: str | None
    protein_id: str | None
    translation: str | None
    gene_kind: str
    gene_functions: list[GeneFunction]
    ec_numbers: list[str]
    db_xrefs: list[str]
    notes: list[str]
    inference: list[str]
    region_numbers: list[int]
    candidate_cluster_numbers: list[int]
    protocluster_numbers: list[int]
    proto_core_numbers: list[int]
    qualifiers: Qualifiers


@dataclass(slots=True)
class CollectionFeature:
    feature_type: str
    number: int | None
    location: Location
    products: list[str]
    references: list[int]
    kind: str | None
    category: str | None
    rules: list[str]
    smiles: list[str]
    polymer: list[str]
    core_location: str | None
    cutoff: int | None
    neighbourhood: int | None
    creating_tool: str | None
    contig_edge: bool | None
    qualifiers: Qualifiers


@dataclass(slots=True)
class Domain:
    domain_id: str | None
    name: str | None
    subtypes: list[str]
    tool: str | None
    locus_tag: str | None
    score: float | None
    evalue: float | None
    protein_start: int | None
    protein_end: int | None
    specificity: list[str]
    location: Location
    qualifiers: Qualifiers

    @property
    def is_nrps_pks(self) -> bool:
        return self.tool == "nrps_pks_domains"


@dataclass(slots=True)
class Module:
    location: Location
    domain_ids: list[str]
    locus_tags: list[str]
    module_type: str | None
    complete: bool
    starter: bool
    final: bool
    iterative: bool
    monomer_pairings: list[str]
    multi_cds: bool
    missing_domain_ids: list[str] = field(default_factory=list)
    qualifiers: Qualifiers = field(default_factory=dict)


@dataclass(slots=True)
class Motif:
    location: Location
    label: str | None
    locus_tag: str | None
    tool: str | None
    score: float | None
    evalue: float | None
    prepeptide: str | None
    core_sequence: str | None
    qualifiers: Qualifiers


@dataclass(slots=True)
class PfamHit:
    location: Location
    accession: str | None
    description: str | None
    locus_tag: str | None
    tool: str | None
    score: float | None
    evalue: float | None
    protein_start: int | None
    protein_end: int | None
    qualifiers: Qualifiers

    def dedup_key(self, record_id: str) -> tuple[Any, ...]:
        normalized = self.accession.split(".", 1)[0] if self.accession else None
        return (
            record_id,
            self.locus_tag,
            self.location.start,
            self.location.end,
            self.protein_start,
            self.protein_end,
            normalized,
        )


@dataclass(slots=True, frozen=True)
class AntiSmashProvenance:
    version: str | None = None
    run_date: str | None = None
    pfam_version: str | None = None
    detection_rule_set_version: str | None = None
    database_versions: dict[str, str] = field(default_factory=dict)
    raw_fields: dict[str, tuple[str, ...]] = field(default_factory=dict)


@dataclass(slots=True)
class Record:
    record_id: str
    name: str
    description: str
    length: int
    molecule_type: str | None
    topology: str | None
    source_path: Path
    source_sha256: str
    antismash_version: str | None
    organism: str | None
    taxonomy: list[str]
    antismash_provenance: AntiSmashProvenance = field(default_factory=AntiSmashProvenance)
    regions: list[CollectionFeature] = field(default_factory=list)
    candidate_clusters: list[CollectionFeature] = field(default_factory=list)
    protoclusters: list[CollectionFeature] = field(default_factory=list)
    proto_cores: list[CollectionFeature] = field(default_factory=list)
    genes: list[Gene] = field(default_factory=list)
    domains: list[Domain] = field(default_factory=list)
    modules: list[Module] = field(default_factory=list)
    motifs: list[Motif] = field(default_factory=list)
    pfam_hits: list[PfamHit] = field(default_factory=list)
    raw_features: list[RawFeature] = field(default_factory=list)
    diagnostics: list[Diagnostic] = field(default_factory=list)
    clusterblast_results: list[ClusterBlastResult] = field(default_factory=list)

    @property
    def nrps_pks_domains(self) -> list[Domain]:
        return [domain for domain in self.domains if domain.is_nrps_pks]

    @property
    def deduplicated_pfam_hits(self) -> list[PfamHit]:
        seen: set[tuple[Any, ...]] = set()
        unique: list[PfamHit] = []
        for hit in self.pfam_hits:
            key = hit.dedup_key(self.record_id)
            if key not in seen:
                seen.add(key)
                unique.append(hit)
        return unique
