---
name: antismash-review
description: Parse, enrich, review, and compare antiSMASH GenBank outputs and result directories with strict Biopython semantics, BGC hierarchy, CDS/domain/motif models, Pfam hits, boundary diagnostics, ClusterBlast sidecar enrichment, and Markdown, JSON, TSV, entity-level TSVs, or comparative review exports.
---

# antiSMASH GenBank review

Use the packaged `inspect` and `compare` commands for analysis and the Python modules for maintenance. Preserve source evidence and qualify biological interpretations conservatively.

## Run the workflow

1. Install the checkout when needed:

   ```bash
   python -m pip install -e .
   ```

2. Inspect a GenBank file or result directory:

   ```bash
   python -m antismash_review inspect region001.gbk
   python -m antismash_review inspect results/ --format json --output review.json
   python -m antismash_review inspect results/ --format gene-tsv --output genes.tsv
   python -m antismash_review inspect results/ --format domain-tsv --output domains.tsv
   python -m antismash_review inspect results/ --format clusterblast-tsv --output clusterblast.tsv
   python -m antismash_review inspect results/ --recursive --format tsv --output review.tsv
   ```

3. Compare two antiSMASH runs or files:

   ```bash
   python -m antismash_review compare run1/ run2/ --format markdown
   python -m antismash_review compare run1/ run2/ --format json --output comparison.json
   python -m antismash_review compare ref.gb alt.gb --match-by coordinate_overlap --assume-shared-coordinate-system
   ```

   Supported record matching modes: `record_id` (default), `record_region`, `single_record`, and `coordinate_overlap`. Coordinate mode requires `--assume-shared-coordinate-system`; `--min-reciprocal-overlap` optionally changes the default 0.80 threshold. Use coordinate matching only when coordinate correspondence is independently established (e.g. re-annotations of identical contigs), not between arbitrary isolates or rebased region files.

4. Use `--lenient` only to retain records when a recognized feature or sidecar cannot be adapted. It does not repair malformed GenBank or switch parsers. Report every emitted diagnostic.

   5. Read [references/semantic-contract.md](references/semantic-contract.md) before interpreting hierarchy, domain/module counts, motifs, boundary status, Pfam totals, ClusterBlast results, or coordinate matching.

6. Build a reproducible cohort matrix from immediate member directories or an explicit TSV manifest:

   ```bash
   python -m antismash_review cohort strains/ --format product-matrix-tsv
   python -m antismash_review cohort strains/ --format product-matrix-tsv --value count
   python -m antismash_review cohort --manifest samples.tsv --format domain-matrix-tsv --cluster-by domain-jaccard --tree-output domains.nwk
   ```

   Directory mode sorts member names. Manifest mode preserves listed order and accepts `sample<TAB>path` rows. Exactly one root or manifest is required. Invalid members fail with their name and path by default; `--skip-invalid-members` is explicit and reported in JSON. Cohort products count regions, domains count adapted `aSDomain` features, and product-class presence is not a homologous BGC-family call.

## Use the Python API

Import the supported library surface from `antismash_review`:

```python
from pathlib import Path

from antismash_review import (
    GenBankParseError,
    dumps_records,
    parse_genbank,
    assess_architecture,
    predict_assembly_lines,
    review_record,
)

try:
    records = parse_genbank(Path("result.gbk"))
except GenBankParseError as exc:
    raise RuntimeError(f"Could not review antiSMASH GenBank: {exc}") from exc

diagnostics = [item for record in records for item in review_record(record)]
json_text = dumps_records(records)
assembly_lines = [prediction for record in records for prediction in predict_assembly_lines(record)]
architecture = [assessment for record in records for assessment in assess_architecture(record)]
```

The stable top-level names are listed in `antismash_review.__all__`. `parse_genbank()` accepts one GenBank file, preserves every record, and raises `GenBankParseError` for unreadable, malformed, or empty input. Its default strict mode raises when a recognized feature cannot be adapted; `lenient=True` retains the record and emits an explicit diagnostic for that feature. It does not discover result directories or attach ClusterBlast sidecars; use the CLI for that enriched workflow.

The supported renderers are `dumps_records()`, `render_records()`, `render_tsv()`, `render_gene_tsv()`, `render_domain_tsv()`, `render_clusterblast_tsv()`, and `render_provenance_json()`. Assembly-line predictions are available through `predict_assembly_lines()` and their dedicated JSON/TSV/Markdown exporters. `build_cohort()` provides the typed cohort-loading API; matrix and clustering presentation is exposed through the CLI/exporter modules. The installed distribution includes `py.typed` so downstream type checkers use the inline annotations.

## Select inputs deliberately

- Accept `.gbk`, `.gb`, and `.gbff` files.
- For a directory, prefer aggregate GenBank files when present; otherwise inspect region files. Discovery is deterministic.
- ClusterBlast sidecars (`clusterblast/`, `knownclusterblast/`, `subclusterblast/`, and native JSON sidecars) are enriched when GenBank input is provided.
- Do not pass native antiSMASH JSON alone as primary input. JSON alone cannot construct the review record model.
- Add `--recursive` only when nested directories belong to the intended input set.

## Interpret results conservatively

- Keep every GenBank record distinct.
- Treat absent `gene_kind` as `unclassified`, not `other`.
- Count an `aSDomain` as NRPS/PKS architecture only when `aSTool=nrps_pks_domains`; retain all other antiSMASH domains separately.
- Keep `aSModule` entities separate and inspect missing or duplicate domain-reference diagnostics.
- Do not apply a universal motif e-value threshold. Interpret motif evidence using the producing tool and motif family.
- Treat domain specificity as source evidence: repeated values, including KR activity or
  stereochemistry strings, are retained but are not independently validated predictions.
- Distinguish boundary-limited context from demonstrated core or CDS truncation.
- Do not infer TFBS validity, LysR identity, HexS orthology, substrate confirmation, frameshifts, or catalytic loss from parser output alone.
- Report both raw and deduplicated Pfam counts when comparing annotations.

## Choose an export

- `markdown`: Concise human review with ClusterBlast tables and review diagnostics.
- `json`: Versioned, structured review envelope (`schema_version: 0.3.0`; includes `antismash_provenance`).
- `tsv`: Compact one-row-per-record summary.
- `gene-tsv`: Entity-level export with one row per gene and zero-based coordinates `[start, end)`.
- `domain-tsv`: Entity-level export with one row per domain and zero-based coordinates.
- `clusterblast-tsv`: Entity-level export with ClusterBlast / KnownClusterBlast / SubClusterBlast hit rankings and pairings.
- `assemblyline-tsv`: Local NRPS/PKS module and monomer calls plus separate Phase 2 core-mass candidate columns; unresolved chemistry remains null (mass candidates generated with v0.1.0 must be regenerated with v0.2.0+).
- `assemblyline-json`: Versioned derived assembly-line evidence and conservative core-mass candidates (`schema_version: 0.3.0`, `parser_version: 0.3.0`).
- `assemblyline-markdown`: Human-readable rendering of the same evidence and caveats.
- `architecture-json`: Scoped domain-slot assessments (`schema_version: 0.1.0`); product classes beyond T1PKS, transAT-PKS, and NRPS are intentionally `not_applicable`; domain evidence is preserved and exported for downstream interpretation.
- `gff3`: Genome-browser GFF3 track with one-based inclusive coordinates and localized findings.
- `bed`: BED6 track retaining zero-based half-open coordinates.
- `provenance-json`: Deduplicated source/run manifest with raw antiSMASH metadata (`schema_version: 0.1.0`).
- `provenance-tsv`: One row per source/hash manifest entry.
- Comparative review (`compare`): Markdown or JSON (`schema_version: 0.2.0`) comparing feature counts, product deltas, diagnostic deltas, intergenic gaps, and provenance deltas.
- `cohort product-matrix-tsv`: Deterministic member-by-product presence/count matrix; products are region-level labels.
- `cohort domain-matrix-tsv`: Deterministic member-by-aSDomain presence/count matrix; Pfam hits are not mixed into this matrix.
- `cohort json`: Versioned cohort envelope (`schema_version: 0.1.0`) with member inputs, source hashes, provenance summaries, normalized/raw column labels, matrices, and explicit skipped members.
- `--cluster-by domain-jaccard`: Optional binary-domain Jaccard distance matrix, deterministic average-linkage leaf order, and Newick tree; it never runs unless explicitly requested.

## Maintain the codebase

Run all checks after changing parsing, models, discovery, exporters, CLI behavior, or this skill.

PowerShell:

```powershell
python -m ruff check .
python -m ruff format --check antismash_review tests
python -m mypy antismash_review
$pytestBase = Join-Path $env:TEMP ("antismash-review-pytest-" + [guid]::NewGuid())
python -m pytest -p no:cacheprovider --basetemp=$pytestBase -q
python -m pytest -p no:cacheprovider --basetemp=$pytestBase `
  --cov=antismash_review --cov-report=term-missing -q
python -m antismash_review --help
```

POSIX shells:

```bash
python -m ruff check .
python -m ruff format --check antismash_review tests
python -m mypy antismash_review
pytestBase="$(mktemp -d -t antismash-review-pytest-XXXXXX)"
python -m pytest -p no:cacheprovider --basetemp="$pytestBase" -q
python -m pytest -p no:cacheprovider --basetemp="$pytestBase" \
  --cov=antismash_review --cov-report=term-missing -q
python -m antismash_review --help
```

Keep local biological integration files private unless redistribution permission is explicit. Tests that use those files must skip cleanly when they are absent.
