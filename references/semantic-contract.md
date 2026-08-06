# Semantic contract

## Parsing and provenance

- Parse the original file with Biopython. Do not pre-reflow text or silently fall back to another parser.
- Calculate a source SHA-256 and retain the resolved source path on each record.
- Return one model object per GenBank record; never flatten multi-record files.
- Retain every feature and its list-valued qualifiers in `raw_features`, including unknown feature types.
- In strict mode, fail when a recognized feature cannot be adapted. In lenient mode, retain a `feature_adapter_failed` diagnostic and continue. Parser-level failures remain fatal.

## Feature semantics

- Preserve `region`, `cand_cluster`, `protocluster`, and `proto_core` as separate collection layers.
- Preserve repeated products, rules, references, SMILES, and polymer values. Normalize whitespace only for machine-valued fields whose tokens can wrap.
- Derive gene collection memberships by coordinate overlap while preserving compound and fuzzy location parts.
- Store missing `gene_kind` as `unclassified`; reserve `other` for an explicit source qualifier.
- Retain all `aSDomain` features in `domains`. Derive `nrps_pks_domains` only from domains whose tool is exactly `nrps_pks_domains`.
- Parse plural `domain_subtypes` and legacy singular `domain_subtype`.
- Keep `aSModule` objects separate from domains. Resolve their domain IDs and diagnose missing or duplicated references.
- Allow `CDS_motif` records without a motif label or e-value, including RiPP prepeptide annotations.
- Retain raw `PFAM_domain` hits and expose a stable deduplicated view keyed by record, locus, nucleotide/protein coordinates, and version-normalized accession.

## Review boundaries

- `region.contig_edge=True` means the extracted BGC context reaches a record boundary. It does not prove that a biosynthetic core or CDS is truncated.
- Diagnose a proto-core at a record edge separately.
- Diagnose a CDS edge truncation only when its location is fuzzy/partial and reaches the edge.
- Do not emit a generic motif-confidence warning based on one e-value cutoff.

## Discovery and exports

- Classify recognized files deterministically as native JSON, aggregate GenBank, or region GenBank.
- During inspection, consume aggregate GenBank files in preference to region files to avoid duplicate representations.
- Native antiSMASH JSON discovery exists only for manifest classification; parsing and source merging require future version-specific adapters.
- Markdown and the top-level JSON diagnostics include evidence-scoped review diagnostics.
- The JSON envelope reports `schema_name`, `schema_version`, and `parser_version`; it is pre-1.0 and may evolve.
- TSV is a compact one-row-per-record summary and is not an entity-level export.
