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

4. Use `--lenient` only to retain records when a recognized feature or sidecar cannot be adapted. It does not repair malformed GenBank or switch parsers. Report every emitted diagnostic.

5. Read [references/semantic-contract.md](references/semantic-contract.md) before interpreting hierarchy, domain/module counts, motifs, boundary status, Pfam totals, ClusterBlast results, or coordinate matching.

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
- Distinguish boundary-limited context from demonstrated core or CDS truncation.
- Do not infer TFBS validity, LysR identity, HexS orthology, substrate confirmation, frameshifts, or catalytic loss from parser output alone.
- Report both raw and deduplicated Pfam counts when comparing annotations.

## Choose an export

- `markdown`: Concise human review with ClusterBlast tables and review diagnostics.
- `json`: Versioned, structured review envelope (`schema_version: 0.2.0`).
- `tsv`: Compact one-row-per-record summary.
- `gene-tsv`: Entity-level export with one row per gene and zero-based coordinates `[start, end)`.
- `domain-tsv`: Entity-level export with one row per domain and zero-based coordinates.
- `clusterblast-tsv`: Entity-level export with ClusterBlast / KnownClusterBlast / SubClusterBlast hit rankings and pairings.
- Comparative review (`compare`): Markdown or JSON (`schema_version: 0.1.0`) comparing feature counts, product deltas, diagnostic deltas, and intergenic gaps.

## Maintain the codebase

Run all checks after changing parsing, models, discovery, exporters, CLI behavior, or this skill:

```bash
python -m ruff check .
python -m ruff format --check antismash_review tests
python -m mypy antismash_review
python -m pytest -q
python -m pytest --cov=antismash_review --cov-report=term-missing
python -m antismash_review --help
```

Keep local biological integration files private unless redistribution permission is explicit. Tests that use those files must skip cleanly when they are absent.
