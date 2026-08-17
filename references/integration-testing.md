# Optional private integration testing

Set `ANTISMASH_REVIEW_INTEGRATION_MANIFEST` to a TSV manifest when private,
redistribution-restricted antiSMASH result directories are available. The suite skips
cleanly when the variable is absent and fails on a configured but missing input or
fingerprint.

The manifest header is:

```text
name	path	fingerprint	recursive	lenient
```

`name`, `path`, and `fingerprint` are required. Relative paths resolve from the manifest
directory. `recursive` and `lenient` are optional booleans accepting `true`/`false`,
`yes`/`no`, or `1`/`0`.

Each fingerprint is exact JSON with `schema_version: "0.1.0"`. It freezes deterministic
record IDs, antiSMASH versions, sidecar basenames/types and schema provenance,
region/product/domain/module totals, assembly-line multiplicity, pairing status, and
integrity flags, mass nullability, diagnostic-code counts, and selected
KnownClusterBlast accessions. It
does not contain absolute paths, timestamps, or nucleotide/protein sequences.

Run the suite with the same full validation command used for ordinary tests. If a
fingerprint changes, inspect the structural JSON diff and update it only when the source
fixture or intended parser semantics changed.
