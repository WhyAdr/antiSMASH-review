---
name: antismash-review
description: Parse and review antiSMASH GenBank files and result directories with strict per-record Biopython semantics, BGC hierarchy, genes, NRPS/PKS domains and modules, motifs, Pfam hits, boundary diagnostics, and Markdown, JSON, or TSV exports. Use for single-region inspection, aggregate GenBank review, deterministic directory discovery, BGC annotation audits, or maintenance of the antismash-review parser and compatibility CLI.
---

# antiSMASH GenBank review

Use the packaged `inspect` command for analysis and the Python modules for maintenance. Preserve source evidence and qualify biological interpretations conservatively.

## Run the workflow

1. Install the checkout when needed:

   ```bash
   python -m pip install -e .
   ```

2. Inspect a GenBank file or result directory:

   ```bash
   python -m antismash_review inspect region001.gbk
   python -m antismash_review inspect results/ --format json --output review.json
   python -m antismash_review inspect results/ --recursive --format tsv --output review.tsv
   ```

3. Use `--lenient` only to retain records when a recognized feature cannot be adapted. It does not repair malformed GenBank or switch parsers. Report every emitted diagnostic.

4. Read [references/semantic-contract.md](references/semantic-contract.md) before interpreting hierarchy, domain/module counts, motifs, boundary status, or Pfam totals.

## Select inputs deliberately

- Accept `.gbk`, `.gb`, and `.gbff` files.
- For a directory, prefer aggregate GenBank files when present; otherwise inspect region files. Discovery is deterministic.
- Add `--recursive` only when nested directories belong to the intended input set.
- Do not pass native antiSMASH JSON as input. JSON is currently an output envelope; version-aware antiSMASH JSON adapters and GenBank/JSON merging are not implemented.

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

- Use Markdown for a concise human review.
- Use JSON for the versioned, structured pre-1.0 review envelope.
- Use TSV for a compact one-row-per-record comparison table.

## Maintain the codebase

Run all checks after changing parsing, models, discovery, exporters, CLI behavior, or this skill:

```bash
python -m ruff check .
python -m mypy antismash_review
python -m pytest -q
python -m pytest --cov=antismash_review --cov-report=term-missing
python -m antismash_review --help
```

Keep local biological integration files private unless redistribution permission is explicit. Tests that use those files must skip cleanly when they are absent.
