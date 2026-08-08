# antiSMASH-review fifth-pass API wrap-up

## Finding

The package already contains a useful importable, programmatic workflow beyond the CLI:

```python
from pathlib import Path

from antismash_review.genbank import parse_genbank
from antismash_review.exporters.json_export import dumps_records

records = parse_genbank(Path("result.gbk"))
json_text = dumps_records(records)
```

The current flow is:

```text
parse_genbank(Path) -> list[Record] -> review/export functions
```

Evidence from the current source:

- `antismash_review.genbank.parse_genbank()` is type-annotated as
  `parse_genbank(path: Path, *, lenient: bool = False) -> list[Record]`.
- `Record` and the related domain, gene, module, motif, Pfam, diagnostic, and ClusterBlast
  models are typed dataclasses in `antismash_review.models`.
- Exporters accept the parsed records, including `dumps_records()`, `render_records()`,
  `render_tsv()`, `render_gene_tsv()`, `render_domain_tsv()`, and
  `render_clusterblast_tsv()`.
- `review_record(record)` provides programmatic review diagnostics.

## Qualification

This is currently an implicit or low-level programmatic API rather than a formally declared
stable public API:

- `antismash_review/__init__.py` currently exposes only `__version__`.
- `from antismash_review import parse_genbank` does not currently work; callers must import
  from the implementation submodules.
- The existing `SKILL.md` documents CLI usage and semantic behavior but does not document
  supported Python imports or programmatic examples.
- There is no `antismash_review/py.typed` marker in the package.

## Meaning of `py.typed`

The `py.typed` marker is a PEP 561 package-data marker. It tells downstream type checkers such
as mypy and Pyright that the installed distribution intentionally ships inline type annotations
that should be checked. The source is already substantially annotated and passes the repository
mypy gate, but the installed package should include the marker for consumers to receive that
typing contract.

Adding the marker requires both:

1. an empty `antismash_review/py.typed` file; and
2. package-data configuration ensuring the marker is included in the wheel.

## Recommended API wrap-up

If the package is intended for reuse as a library, add a small, explicit API layer:

1. Document supported imports and a minimal parse-review-export example.
2. Decide which names are stable public imports.
3. Optionally re-export the selected names from `antismash_review.__init__`, for example
   `parse_genbank`, `GenBankParseError`, `Record`, `review_record`, and selected exporters.
4. Add `__all__` for the chosen public names rather than making every submodule symbol an API
   promise.
5. Add and package `antismash_review/py.typed`.
6. Document strict versus lenient parsing and the public exception behavior.
7. Add import-level API smoke tests and verify the marker is present in the built wheel.

These changes would improve discoverability and downstream typing without changing the existing
GenBank parsing semantics or CLI behavior.
