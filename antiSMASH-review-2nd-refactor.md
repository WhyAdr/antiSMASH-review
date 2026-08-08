# antiSMASH-review second refactor — reviewed plan

This revision is grounded in the current package, its 19-test baseline, and the local antiSMASH 8.0.4 result bundle `SM-ZPG19--NOACCESSION-antismash`. It keeps the parser GenBank-first and adds sidecar enrichment without pretending that a ClusterBlast-only JSON adapter is a complete native antiSMASH JSON adapter.

## Review outcome

The earlier plan had several sound goals, but it was not safe to implement verbatim.

| Earlier proposal | Finding | Refined decision |
|---|---|---|
| “All changes are additive” | Adding fields to JSON and columns to TSV changes public schemas. Accepting JSON enrichment also changes directory behavior. | Version every changed schema and call out compatibility changes explicitly. |
| Inspect native JSON by itself | The CLI still constructs `Record` objects only from GenBank. ClusterBlast JSON cannot supply the rest of the record model. | Native JSON remains unsupported as the primary input. JSON is an optional sidecar beside GenBank. |
| Flat `ClusterBlastHit` list | It loses empty results, per-region grouping, source provenance, and JSON `total_hits`. | Model a `ClusterBlastResult` containing ranked hits. |
| Attach hits by `record_id` alone | Ten region GBKs all use `contig_1`; every result would be attached to every region record. | Match on record identity **and** `region_number`, requiring exactly one target. |
| Prefer all text if any text exists | One text file would suppress JSON fallback for every missing mode and region. | Prefer text per `(record_id, region_number, search_type)` key; fill only missing keys from JSON. |
| Swallow JSON exceptions | Silent loss violates the strict parser contract. | Raise a public enrichment error; `--lenient` may convert it to an explicit diagnostic. |
| `clusterblast_hits` count | In SM-ZPG19, general JSON reports 48,710 database hits but only 500 ranked rows. | Name counts `total_hits` and `ranked_hit_count`; never conflate them. |
| High Pfam duplication diagnostic at 25% | All ten real region GBKs deduplicate to about 50%, so this would emit routine false positives. | Do not add this diagnostic. Continue reporting raw and deduplicated counts. |
| Match different isolates by absolute coordinate overlap | Coordinates from unrelated assemblies can overlap numerically without homology. | Replace `locus_position` with guarded `coordinate_overlap`, requiring an explicit shared-coordinate-system assumption. |
| Compare consecutive `Location.start/end` values | Compound/cross-origin CDSs span misleading min/max bounds; circular wrap gaps are omitted. | Calculate an interval-union summary from location parts and topology. |

## Ground-truth baselines

Current unit baseline:

- `19 passed`.
- Overall coverage is `93%`.
- The existing `.pytest_cache` is not writable in this checkout, so validation must use `-p no:cacheprovider` to avoid a known warning.

SM-ZPG19 aggregate GenBank baseline:

| Metric | Value |
|---|---:|
| antiSMASH version | 8.0.4 |
| Records | 1 (`contig_1`) |
| Length | 5,269,270 bp |
| Regions | 10 |
| CDS features | 4,877 |
| All `aSDomain` features | 178 |
| NRPS/PKS domains | 56 |
| Modules | 10 |
| Motifs | 109 |
| Pfam hits | 6,967 raw / 6,548 deduplicated |

SM-ZPG19 ClusterBlast baseline:

| Search type | Results | Non-empty results | Ranked rows | JSON `total_hits` |
|---|---:|---:|---:|---:|
| ClusterBlast | 10 | 10 | 500 | 48,710 |
| KnownClusterBlast | 10 | 7 | 63 | 63 |
| SubClusterBlast | 10 | 0 | 0 | 0 |

The accession, description, and rank sequences in all 30 text files agree exactly with the corresponding JSON rankings. Text files do not expose every JSON score field, particularly similarity. Treat `total_hits` as an antiSMASH-provided counter whose precise upstream meaning belongs to the supported schema; do not relabel it as the number of displayed rankings.

The local SM-ZPG19 folder is private integration evidence. Keep it untracked and do not copy it into `tests/fixtures` or a distributable artifact.

## Scope and phase gates

Each phase must pass its focused tests and the full quality gate before the next begins. Phase 4 is deliberately gated on the comparison contract; coordinate matching is available only when the user explicitly asserts that both inputs share a coordinate system.

### Phase 1 — Focused coverage and conservative diagnostics

#### Housekeeping

Update [.gitignore](.gitignore) before running or staging tests:

```gitignore
/.pytest_temp/
/SM-ZPG19--NOACCESSION-antismash/
```

The first entry covers the checkout-local test base used below. The second protects this specific private integration bundle from accidental staging; it does not redistribute or generalize the biological input.

#### Test fixtures

Add [tests/conftest.py](tests/conftest.py) with one reusable synthetic-fixture root and one optional private-integration fixture. Avoid a separate pytest fixture for every static filename.

```python
from __future__ import annotations

import os
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture
def sm_zpg19_dir() -> Path:
    configured = os.environ.get("ANTISMASH_SM_ZPG19_DIR")
    candidates = [
        Path(configured).expanduser() if configured else None,
        REPOSITORY_ROOT / "SM-ZPG19--NOACCESSION-antismash",
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_dir():
            return candidate.resolve()
    pytest.skip("private SM-ZPG19 antiSMASH result directory is unavailable")
```

Add two small synthetic GenBank fixtures:

- [tests/fixtures/partial-edge.gb](tests/fixtures/partial-edge.gb): a linear 100 bp record with a fuzzy CDS at `<1..30`.
- [tests/fixtures/cross-origin.gb](tests/fixtures/cross-origin.gb): a circular 100 bp record with a compound CDS spanning the origin.

Add or relocate tests so ownership is clear:

- [tests/test_genbank.py](tests/test_genbank.py): empty GenBank raises `GenBankParseError`; cross-origin location parts and collection membership are retained.
- [tests/test_review.py](tests/test_review.py): partial edge, orphan module, and missing NRPS/PKS architecture diagnostics.
- [tests/test_exporters.py](tests/test_exporters.py): all `contig_edge` tri-state outputs and Markdown/TSV escaping.
- [tests/test_cli.py](tests/test_cli.py): no-command behavior, discovery errors, parse errors, output-write errors, overwrite protection, and JSON-only rejection.
- [tests/test_discovery.py](tests/test_discovery.py): unsupported single files and nonexistent directories.

#### Diagnostics

Add only two consistency diagnostics to [antismash_review/review.py](antismash_review/review.py). They identify mismatches in parsed antiSMASH evidence; they do not make biological-quality claims.

```python
def _extend_consistency_diagnostics(
    record: Record,
    diagnostics: list[Diagnostic],
) -> None:
    gene_locus_tags = {gene.locus_tag for gene in record.genes if gene.locus_tag}
    for module_index, module in enumerate(record.modules):
        orphan_tags = sorted({tag for tag in module.locus_tags if tag not in gene_locus_tags})
        if orphan_tags:
            diagnostics.append(
                Diagnostic(
                    code="orphan_module_locus",
                    severity=Severity.WARNING,
                    message=(
                        f"Module {module_index + 1} references locus tags absent from "
                        f"the CDS set: {', '.join(orphan_tags)}"
                    ),
                    source=str(record.source_path),
                    record_id=record.record_id,
                )
            )

    architecture_products = sorted(
        {
            product
            for region in record.regions
            for product in region.products
            if "nrps" in product.casefold() or "pks" in product.casefold()
        }
    )
    if architecture_products and not record.nrps_pks_domains:
        diagnostics.append(
            Diagnostic(
                code="missing_nrps_pks_architecture",
                severity=Severity.WARNING,
                message=(
                    "Region products imply NRPS/PKS architecture "
                    f"({', '.join(architecture_products)}), but no domains from "
                    "aSTool=nrps_pks_domains were parsed"
                ),
                source=str(record.source_path),
                record_id=record.record_id,
            )
        )
```

Call this helper once at the end of `review_record()`. Do not add `high_pfam_duplication`; raw and deduplicated Pfam totals already expose the underlying representation.

Document every existing and new diagnostic in [references/semantic-contract.md](references/semantic-contract.md), including severity, trigger, and non-claim.

#### Phase 1 acceptance

- Existing outputs remain byte-for-byte unchanged unless one of the two new diagnostics legitimately fires.
- The SM-ZPG19 aggregate record fires neither new diagnostic.
- Region-level Pfam duplication does not create a diagnostic.
- Full suite and static checks pass.

### Phase 2 — Entity-level TSV exports

Add [antismash_review/exporters/entity_tables.py](antismash_review/exporters/entity_tables.py) with `render_gene_tsv()` and `render_domain_tsv()`. Use `csv.writer`; never assemble rows with string concatenation.

Shared serialization rules:

```python
def _list_cell(values: list[str] | list[int]) -> str:
    return "; ".join(str(value) for value in values)


def _optional(value: object | None) -> object:
    return "" if value is None else value


def _boolean(value: bool) -> str:
    return "true" if value else "false"
```

Gene TSV columns, in stable order:

```python
GENE_COLUMNS = (
    "source_path",
    "source_sha256",
    "filename",
    "record_id",
    "locus_tag",
    "gene",
    "product",
    "gene_kind",
    "start",
    "end",
    "strand",
    "partial",
    "cross_origin",
    "ec_numbers",
    "db_xrefs",
    "region_numbers",
    "candidate_cluster_numbers",
    "protocluster_numbers",
    "proto_core_numbers",
)
```

Domain TSV columns, in stable order:

```python
DOMAIN_COLUMNS = (
    "source_path",
    "source_sha256",
    "filename",
    "record_id",
    "domain_id",
    "name",
    "tool",
    "is_nrps_pks",
    "locus_tag",
    "start",
    "end",
    "strand",
    "partial",
    "cross_origin",
    "protein_start",
    "protein_end",
    "score",
    "evalue",
    "subtypes",
    "specificity",
)
```

Coordinates remain the model's zero-based, half-open coordinates. State that explicitly in the semantic contract and TSV documentation; do not silently convert them to GenBank's one-based display coordinates.

Extend `inspect --format` in [antismash_review/cli.py](antismash_review/cli.py):

```python
inspect.add_argument(
    "--format",
    choices=("markdown", "json", "tsv", "gene-tsv", "domain-tsv"),
    default="markdown",
)
```

Tests must assert exact headers, row counts, empty optional values, lowercase booleans, list ordering, tabs/newlines in products, cross-origin fields, and an NRPS domain versus a non-NRPS TIGRFAM domain.

#### Phase 2 acceptance

- Existing `tsv` remains the compact one-row-per-record schema and does not change in this phase.
- Gene and domain TSVs parse with `csv.reader(delimiter="\t")` and have a constant number of cells per row.
- The two new CLI formats work both on stdout and with `--output`.

### Phase 3 — Version-gated ClusterBlast sidecar enrichment

#### Models and schema

Add the following models to [antismash_review/models.py](antismash_review/models.py). Preserve antiSMASH's public `search_type` names instead of inventing `general`/`knowncluster` aliases.

```python
from typing import Literal

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


# Add to Record:
clusterblast_results: list[ClusterBlastResult] = field(default_factory=list)
```

Create [antismash_review/schema.py](antismash_review/schema.py):

```python
RECORD_SCHEMA_NAME = "antismash-review"
RECORD_SCHEMA_VERSION = "0.2.0"
COMPARISON_SCHEMA_NAME = "antismash-review-comparison"
COMPARISON_SCHEMA_VERSION = "0.1.0"
```

Use these constants in JSON exporters. The record schema moves from `0.1.0` to `0.2.0` when `clusterblast_results` is added; this is not a no-op refactor.

#### Discovery

Extend `InputManifest` in [antismash_review/discovery.py](antismash_review/discovery.py):

```python
@dataclass(slots=True, frozen=True)
class InputManifest:
    root: Path
    json_files: tuple[Path, ...]
    aggregate_genbanks: tuple[Path, ...]
    region_genbanks: tuple[Path, ...]
    ignored_files: tuple[Path, ...]
    clusterblast_files: tuple[Path, ...] = ()
    knownclusterblast_files: tuple[Path, ...] = ()
    subclusterblast_files: tuple[Path, ...] = ()
```

For a result-directory input, inspect only the canonical immediate children `clusterblast/*.txt`, `knownclusterblast/*.txt`, and `subclusterblast/*.txt`. Sort paths by a natural `(contig, region_number)` key so `c2` precedes `c10`. A single-file input has no automatically discovered sidecars.

Do not claim that the current flat `InputManifest` safely groups multiple recursively nested antiSMASH runs. Recursive multi-run grouping is a separate future refactor.

#### Parsers

Create [antismash_review/clusterblast.py](antismash_review/clusterblast.py):

```python
class ClusterBlastParseError(RuntimeError):
    """A ClusterBlast sidecar was recognized but could not be parsed safely."""


def parse_clusterblast_text(
    path: Path,
    *,
    search_type: ClusterBlastSearchType,
) -> ClusterBlastResult:
    """Parse one antiSMASH 8 ClusterBlast text result."""
    ...


def parse_clusterblast_json(path: Path) -> list[ClusterBlastResult]:
    """Parse only ClusterBlast modules from a native antiSMASH JSON sidecar."""
    ...
```

Text parser contract:

1. Parse `ClusterBlast scores for <record_id>`; all three modes use that same header.
2. Extract `region_number` from `_cN.txt`; a non-matching filename is an error, never region `0`.
3. Accept trailing whitespace on `Significant hits:` and `Details:` markers.
4. Require consecutive one-based ranks and make the details accession agree with the significant-hit accession.
5. Parse exactly six tab-separated BLAST columns after the BLAST-table header.
6. Return an empty `rankings` list for a valid empty result.
7. Set text `total_hits=None`; the text file exposes ranked results, not the JSON `total_hits` counter.
8. Retain path and SHA-256 provenance.

JSON parser contract:

1. Read only `records[*].modules["antismash.modules.clusterblast"]`.
2. A valid JSON document without that module returns an empty list.
3. Validate the observed module schema `2` and result schema `5`. Unsupported schemas raise `ClusterBlastParseError` with the path and observed values.
4. Validate module `record_id` against the containing JSON record ID.
5. Preserve `region_number`, `total_hits`, all ranking score fields, pairings, source path, SHA-256, and schema versions.
6. Extract the query locus tag with `query_string.split("|", 5)[4]`; reject malformed query strings.
7. Do not catch `KeyError` or `ValueError` with `pass`.

#### Source selection and attachment

Source precedence is per result key, not global:

```python
def _result_key(
    result: ClusterBlastResult,
) -> tuple[str, int, ClusterBlastSearchType]:
    return (result.record_id, result.region_number, result.search_type)


def merge_clusterblast_results(
    text_results: list[ClusterBlastResult],
    json_results: list[ClusterBlastResult],
) -> list[ClusterBlastResult]:
    selected: dict[
        tuple[str, int, ClusterBlastSearchType],
        ClusterBlastResult,
    ] = {}
    for result in json_results:
        key = _result_key(result)
        if key in selected:
            raise ClusterBlastParseError(f"duplicate JSON ClusterBlast result: {key}")
        selected[key] = result
    for result in text_results:
        key = _result_key(result)
        if key in selected and selected[key].source_format == "text":
            raise ClusterBlastParseError(f"duplicate text ClusterBlast result: {key}")
        selected[key] = result
    return [selected[key] for key in sorted(selected)]
```

Attach each result to exactly one parsed record:

```python
def attach_clusterblast_results(
    records: list[Record],
    results: list[ClusterBlastResult],
) -> None:
    for result in results:
        candidates = [
            record
            for record in records
            if result.record_id in {record.record_id, record.name}
            and any(region.number == result.region_number for region in record.regions)
        ]
        if len(candidates) != 1:
            raise ClusterBlastParseError(
                "expected one GenBank target for "
                f"{_result_key(result)}, found {len(candidates)}"
            )
        candidates[0].clusterblast_results.append(result)
```

This works for the preferred aggregate GBK and for a directory containing only region GBKs. In strict mode any recognized malformed or unattached sidecar is fatal. In `--lenient` mode, convert the failure into a `clusterblast_parse_failed` diagnostic; never discard it silently.

JSON alone remains an error in the CLI:

```text
native antiSMASH JSON cannot yet provide the review record model; provide GenBank,
optionally in a result directory with JSON used as ClusterBlast enrichment
```

#### Exports

Markdown must group by region and search type, not merely by search type:

```markdown
### ClusterBlast

#### Region 1 — KnownClusterBlast

- Source: text (`knownclusterblast/contig_1_c1.txt`)
- Ranked hits: 32

| Rank | Accession | Description | Proteins hit | Score | Similarity |
|---:|---|---|---:|---:|---:|
| 1 | BGC0002414.3 | trichrysobactin/... | 6 | 3338.0 |  |
```

Escape pipes, tabs, and newlines in Markdown table cells. Show the first five rankings per result and state when rows are omitted. Display JSON `total_hits` separately when present.

Add `clusterblast-tsv` as a third entity export. Its stable columns are:

```python
CLUSTERBLAST_COLUMNS = (
    "record_id",
    "region_number",
    "search_type",
    "result_source",
    "result_source_path",
    "total_hits",
    "ranked_hit_count",
    "rank",
    "accession",
    "description",
    "cluster_type",
    "num_hits",
    "blast_score",
    "similarity",
    "pairing_count",
)
```

Emit one blank-rank row for a valid empty result so that negative SubClusterBlast evidence is retained. Do not add an ambiguous `clusterblast_hits` column to the compact record TSV. If compact columns are later desired, name them explicitly (`clusterblast_ranked_rows`, `clusterblast_total_hits`) and bump a documented TSV schema.

#### Tests

Synthetic tests in [tests/test_clusterblast.py](tests/test_clusterblast.py):

- Full text result, empty result, malformed filename, non-consecutive ranks, mismatched details accession, malformed BLAST row, and numeric conversion failure.
- JSON module absent, supported schemas, unsupported module/result schemas, malformed query ID, multiple records, all three search types, empty results, and score-field preservation.
- Text preference for one key while JSON fills another missing key.
- Aggregate attachment and region-only attachment.
- Ambiguous/unattached results fail explicitly.
- Markdown escaping and `clusterblast-tsv` empty-result rows.

Private integration tests in [tests/test_sm_zpg19_integration.py](tests/test_sm_zpg19_integration.py):

- Assert the GenBank baseline table above.
- Assert 10/10/10 text results with 500/63/0 ranked rows.
- Assert JSON 10/10/10 results, 48,710/63/0 total hits, and 500/63/0 ranked rows.
- Assert text and JSON rank/accession/description parity for all 30 result keys.
- Assert a directory with aggregate GBK plus JSON but no text uses JSON fallback.
- Assert JSON by itself still returns CLI status 2.

These integration tests must skip cleanly when the private directory is absent.

#### Phase 3 acceptance

- `inspect SM-ZPG19--NOACCESSION-antismash` attaches exactly 30 results to the one aggregate record without duplicates.
- A region-only copy attaches one result per `(record_id, region_number, search_type)` to the correct region record.
- Text is selected for available keys; JSON fills only missing keys.
- JSON output reports record schema `0.2.0`.
- No supported sidecar parse failure is silent.

### Phase 4 — Comparative review, with explicit record matching semantics

This phase compares parsed evidence; its record matching modes do not independently infer homology.

Supported record matching modes:

1. `record_id` (default): exact ID match; IDs must be unique on each side.
2. `record_region`: key by `(record_id, sole_region_number)`; valid only when every record contains exactly one numbered region.
3. `single_record`: require exactly one record on each side and pair them explicitly even if IDs differ.
4. `coordinate_overlap`: match records by reciprocal overlap of their parsed feature spans. This mode requires `--assume-shared-coordinate-system` and rejects ambiguous or weak matches.

`coordinate_overlap` is intended for the same assembly, re-annotations of the same contigs, or other inputs whose coordinate correspondence is already known. It is not intended to discover homology between arbitrary isolates. In particular, antiSMASH region GBKs commonly rebase each extracted region to coordinate zero; those files do **not** share their original chromosome coordinate system merely because they came from the same run.

The matching algorithm must:

1. Build each record's feature span from the union of collection and CDS location parts.
2. Calculate reciprocal overlap: overlap divided by the left span and by the right span.
3. Require both reciprocal overlaps to meet a documented threshold (default `0.80`).
4. Select one-to-one matches only; reject ties and any case where two left records select the same right record.
5. Report the overlap bases and both reciprocal fractions in the comparison output.

The assumption flag acknowledges coordinate correspondence; it does not convert coordinate overlap into evidence of sequence homology. If cross-isolate BGC homology is needed, add sequence/synteny-aware region matching later or require an explicit mapping file.

#### Comparison models

Create [antismash_review/compare.py](antismash_review/compare.py):

```python
@dataclass(slots=True, frozen=True)
class IntergenicSummary:
    gap_count: int
    total_bp: int
    mean_bp: float | None
    median_bp: float | None
    max_bp: int | None
    circular_wrap_included: bool


@dataclass(slots=True, frozen=True)
class DiagnosticFingerprint:
    code: str
    severity: str
    message: str
    feature_index: int | None


@dataclass(slots=True, frozen=True)
class CoordinateMatchEvidence:
    overlap_bp: int
    left_span_bp: int
    right_span_bp: int
    left_overlap_fraction: float
    right_overlap_fraction: float


@dataclass(slots=True)
class RecordComparison:
    left_record_id: str
    right_record_id: str
    match_key: str
    left_region_count: int
    right_region_count: int
    left_gene_count: int
    right_gene_count: int
    left_domain_count: int
    right_domain_count: int
    left_nrps_pks_count: int
    right_nrps_pks_count: int
    gained_products: list[str]
    lost_products: list[str]
    new_diagnostics: list[DiagnosticFingerprint]
    resolved_diagnostics: list[DiagnosticFingerprint]
    left_intergenic: IntergenicSummary
    right_intergenic: IntergenicSummary
    coordinate_evidence: CoordinateMatchEvidence | None = None


@dataclass(slots=True)
class ComparisonResult:
    left_input: Path
    right_input: Path
    match_method: str
    shared_coordinate_system_assumed: bool
    min_reciprocal_overlap: float | None
    matched: list[RecordComparison]
    unmatched_left: list[str]
    unmatched_right: list[str]
```

Use interval parts, not each gene's min/max span, for intergenic summaries:

```python
def intergenic_summary(record: Record) -> IntergenicSummary:
    intervals = sorted(
        (part.start, part.end)
        for gene in record.genes
        for part in gene.location.parts
    )
    if not intervals:
        return IntergenicSummary(0, 0, None, None, None, False)

    merged: list[list[int]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    gaps = [
        right[0] - left[1]
        for left, right in zip(merged, merged[1:], strict=False)
        if right[0] > left[1]
    ]
    circular = record.topology is not None and record.topology.casefold() == "circular"
    if circular:
        wrap_gap = (record.length - merged[-1][1]) + merged[0][0]
        if wrap_gap > 0:
            gaps.append(wrap_gap)

    return IntergenicSummary(
        gap_count=len(gaps),
        total_bp=sum(gaps),
        mean_bp=statistics.fmean(gaps) if gaps else None,
        median_bp=statistics.median(gaps) if gaps else None,
        max_bp=max(gaps) if gaps else None,
        circular_wrap_included=circular,
    )
```

For linear records, this intentionally measures internal gaps only. For circular records, include the origin-spanning gap. It is a descriptive profile; do not align individual gaps across isolates.

Compare products as `Counter` multisets with deterministic sorted expansion. Compare diagnostics as counters of `(code, severity, message, feature_index)`; excluding source path and record ID allows an explicitly paired record with a renamed ID to resolve correctly while retaining repeated diagnostics.

#### CLI and exporters

Add:

```python
compare_cmd = subparsers.add_parser("compare", help="compare two review inputs")
compare_cmd.add_argument("left", type=Path)
compare_cmd.add_argument("right", type=Path)
compare_cmd.add_argument(
    "--match-by",
    choices=("record_id", "record_region", "single_record", "coordinate_overlap"),
    default="record_id",
)
compare_cmd.add_argument(
    "--assume-shared-coordinate-system",
    action="store_true",
    help=(
        "confirm that left and right coordinates are directly comparable; "
        "required with --match-by coordinate_overlap"
    ),
)
compare_cmd.add_argument(
    "--min-reciprocal-overlap",
    type=float,
    default=0.80,
    metavar="FRACTION",
    help="minimum overlap fraction on both records in coordinate mode (default: 0.80)",
)
compare_cmd.add_argument("--format", choices=("markdown", "json"), default="markdown")
compare_cmd.add_argument("--output", type=Path)
compare_cmd.add_argument("--lenient", action="store_true")
compare_cmd.add_argument("--recursive", action="store_true")
```

Validate the option combination before loading inputs:

```python
if args.match_by == "coordinate_overlap":
    if not args.assume_shared_coordinate_system:
        parser.error(
            "--match-by coordinate_overlap requires "
            "--assume-shared-coordinate-system"
        )
    if not 0 < args.min_reciprocal_overlap <= 1:
        parser.error("--min-reciprocal-overlap must be in the interval (0, 1]")
elif args.assume_shared_coordinate_system:
    parser.error(
        "--assume-shared-coordinate-system is only valid with "
        "--match-by coordinate_overlap"
    )
```

Refactor common input loading into one helper used by `inspect` and `compare`; this prevents the two commands from drifting on aggregate preference, lenient parsing, or ClusterBlast enrichment. Overwrite protection must check every GenBank and JSON input on both sides.

Markdown shows left, right, and delta for counts; product multiplicity; diagnostic fingerprints; and descriptive intergenic summaries. It must label `single_record` as an explicit user-requested pairing. For `coordinate_overlap`, it must print the asserted assumption, threshold, overlap bases, and both reciprocal fractions. JSON uses `COMPARISON_SCHEMA_NAME` and `COMPARISON_SCHEMA_VERSION`.

Tests in [tests/test_compare.py](tests/test_compare.py):

- Identical records produce zero deltas.
- Product gains/losses preserve multiplicity and deterministic order.
- Duplicate record IDs are rejected in `record_id` mode.
- `record_region` accepts unique single-region records and rejects aggregate/missing-number records.
- `single_record` pairs different IDs and rejects either side with zero or multiple records.
- `coordinate_overlap` is rejected unless `--assume-shared-coordinate-system` is present.
- Coordinate mode accepts a unique match above the reciprocal threshold and records its evidence.
- Coordinate mode rejects below-threshold overlaps, equal-score ties, and non-one-to-one assignments.
- `--min-reciprocal-overlap` rejects values outside `(0, 1]`.
- The assumption flag is rejected with every non-coordinate matching mode.
- Diagnostic comparisons retain repeated distinct messages.
- Linear interval-union gaps ignore overlapping CDSs.
- Circular gaps include the wrap interval and handle a cross-origin CDS.
- Unmatched IDs are deterministic.
- Comparison JSON is versioned and serializable.
- CLI overwrite protection covers both inputs.

#### Phase 4 acceptance

- Documentation and output call these “record matching modes,” not “identity modes.”
- Coordinate mode requires the explicit assumption flag and never claims cross-isolate homology.
- Duplicate or ambiguous matches fail with status 2 and a clear message.
- Intergenic summaries are correct for the synthetic linear and circular fixtures.
- Same-input comparison has zero count/product/diagnostic deltas.

### Phase 5 — Documentation, compatibility, and release audit

Update [SKILL.md](SKILL.md) and [references/semantic-contract.md](references/semantic-contract.md) only after behavior is implemented and tested.

Document:

- GenBank remains the authoritative record source.
- Native JSON alone is unsupported; supported ClusterBlast schemas may enrich a GenBank result directory.
- Text-over-JSON precedence is per result key.
- Empty ClusterBlast results and source provenance are retained.
- Entity TSV coordinates are zero-based and half-open.
- Compact record TSV versus gene/domain/ClusterBlast entity TSV purposes.
- The four explicit record matching modes, the coordinate assumption flag, and the reciprocal-overlap threshold.
- Why rebased region GBKs and arbitrary isolates must not use coordinate matching without independently established correspondence.
- Every diagnostic trigger and limitation.
- `misc_feature` entries are retained raw; do not infer that they are ClusterBlast evidence.

Review [agents/openai.yaml](agents/openai.yaml) for command/description drift, but keep it unless validation demonstrates a concrete problem.

## File-change summary

| Phase | File | Action |
|---|---|---|
| 1 | `.gitignore` | Ignore the local pytest base and private SM-ZPG19 bundle |
| 1 | `tests/conftest.py` | Add shared and optional private fixtures |
| 1 | `tests/fixtures/partial-edge.gb` | Add synthetic fixture |
| 1 | `tests/fixtures/cross-origin.gb` | Add synthetic fixture |
| 1 | `tests/test_genbank.py` | Add parser edge tests |
| 1 | `tests/test_review.py` | Add focused diagnostic tests |
| 1 | `tests/test_exporters.py` | Add exporter edge tests |
| 1 | `tests/test_cli.py`, `tests/test_discovery.py` | Close error-path gaps |
| 1 | `antismash_review/review.py` | Add two consistency diagnostics |
| 1 | `references/semantic-contract.md` | Document diagnostics |
| 2 | `antismash_review/exporters/entity_tables.py` | Add gene/domain TSVs |
| 2 | `antismash_review/cli.py` | Route entity formats |
| 2 | `tests/test_entity_tables.py` | Verify stable schemas and escaping |
| 3 | `antismash_review/models.py` | Add result/hit/pairing models and record field |
| 3 | `antismash_review/schema.py` | Centralize record/comparison schema versions |
| 3 | `antismash_review/clusterblast.py` | Add strict text and JSON sidecar parsers |
| 3 | `antismash_review/discovery.py` | Discover canonical sidecar directories |
| 3 | `antismash_review/cli.py` | Merge and attach sidecar results |
| 3 | `antismash_review/exporters/markdown.py` | Add per-region ClusterBlast sections |
| 3 | `antismash_review/exporters/entity_tables.py` | Add ClusterBlast TSV |
| 3 | `antismash_review/exporters/json_export.py` | Emit schema `0.2.0` |
| 3 | `tests/test_clusterblast.py` | Add unit/integration-path tests |
| 3 | `tests/test_sm_zpg19_integration.py` | Add optional real-data regressions |
| 4 | `antismash_review/compare.py` | Add explicit record-matching comparison engine |
| 4 | `antismash_review/exporters/compare_markdown.py` | Add Markdown comparison export |
| 4 | `antismash_review/exporters/compare_json.py` | Add versioned comparison JSON |
| 4 | `antismash_review/cli.py` | Add `compare` command and shared loader |
| 4 | `tests/test_compare.py` | Verify matching and intergenic semantics |
| 5 | `SKILL.md`, `references/semantic-contract.md`, `agents/openai.yaml` | Synchronize and audit documentation |

## Quality gate after every phase

```powershell
python -m ruff check .
python -m ruff format --check antismash_review tests
python -m mypy antismash_review
python -m pytest -p no:cacheprovider --basetemp=.pytest_temp -q
python -m pytest -p no:cacheprovider --basetemp=.pytest_temp `
  --cov=antismash_review --cov-report=term-missing -q
git -c safe.directory='D:/W/Skills Claude/antiSMASH-review' diff --check
```

Coverage must not fall below the current 93% baseline without an explicit, reviewed reason.

CLI smoke tests after the relevant phases:

```powershell
python -m antismash_review --help
python -m antismash_review inspect tests/fixtures/semantics.gb
python -m antismash_review inspect tests/fixtures/semantics.gb --format json
python -m antismash_review inspect tests/fixtures/semantics.gb --format tsv
python -m antismash_review inspect tests/fixtures/semantics.gb --format gene-tsv
python -m antismash_review inspect tests/fixtures/semantics.gb --format domain-tsv
python -m antismash_review inspect SM-ZPG19--NOACCESSION-antismash --format clusterblast-tsv
python -m antismash_review compare tests/fixtures/semantics.gb `
  tests/fixtures/semantics.gb --match-by record_id
python -m antismash_review compare tests/fixtures/semantics.gb `
  tests/fixtures/semantics.gb --match-by coordinate_overlap `
  --assume-shared-coordinate-system
python -m antismash_review compare --help
```

Final artifact checks:

- All Markdown files decode as UTF-8 without replacement characters.
- Fenced code blocks are balanced.
- Local Markdown links resolve.
- No trailing whitespace or mojibake is introduced.
- The private SM-ZPG19 directory and generated outputs remain untracked.
- Build/import/entry-point validation and `pip check` pass in an isolated environment.

## Deferred work

- Full native antiSMASH JSON-to-`Record` parsing.
- Recursive grouping of multiple antiSMASH result bundles.
- Homology- or synteny-aware BGC matching across isolates or assembly versions.
- Gap sequence retention and content-level comparison.
- Merging text and JSON fields within the same ClusterBlast result; phase 3 selects one authoritative source per key.
