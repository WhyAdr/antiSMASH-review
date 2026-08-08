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

## Diagnostics reference

| Code | Severity | Trigger | Non-claim |
|---|---|---|---|
| `context_reaches_record_edge` | NOTICE | Any region has `contig_edge=True` | Does not prove a core or CDS is truncated |
| `core_reaches_record_edge` | WARNING | A proto-core location starts at 0 or ends at record length | Does not prove the core gene is incomplete |
| `partial_cds_at_edge` | WARNING | A CDS with a fuzzy/partial location reaches a record boundary | Does not assess the biological completeness of the protein |
| `feature_adapter_failed` | WARNING (lenient) / ERROR (strict) | A recognized feature could not be adapted to its typed model | In strict mode, halts parsing |
| `module_domain_missing` | WARNING | An aSModule references a domain_id absent from the record | Does not assess module completeness |
| `domain_id_duplicated` | WARNING | A domain_id appears more than once across aSDomain features | Does not determine which instance is canonical |
| `orphan_module_locus` | WARNING | An aSModule references locus tags absent from the CDS set | Does not prove the module is misannotated; the CDS may have been filtered |
| `missing_nrps_pks_architecture` | WARNING | Region products contain NRPS/PKS terms but no aSTool=nrps_pks_domains domains were parsed | Does not prove the annotation is wrong; the relevant domains may be in a different record |
| `clusterblast_parse_failed` | WARNING (lenient) | A recognized ClusterBlast text or JSON sidecar could not be parsed or attached safely | Does not repair or silently discard the malformed sidecar; strict mode fails instead |

Region-level Pfam duplication does not create a diagnostic. Raw and deduplicated Pfam totals already expose the underlying representation.

## Discovery and exports

- Classify recognized files deterministically as native JSON, aggregate GenBank, or region GenBank.
- During inspection, consume aggregate GenBank files in preference to region files to avoid duplicate representations.
- Native antiSMASH JSON alone remains unsupported as primary input; when a GenBank input is provided alongside JSON or sidecar text directories, sidecar enrichment parses and attaches ClusterBlast results.
- Markdown and the top-level JSON diagnostics include evidence-scoped review diagnostics.
- The JSON envelope reports `schema_name: "antismash-review"`, `schema_version: "0.2.0"`, and `parser_version`.
- TSV (`--format tsv`) is a compact one-row-per-record summary.
- Entity-level TSV exports (`--format gene-tsv`, `--format domain-tsv`, and `--format clusterblast-tsv`) export one row per entity. All entity TSV coordinates (`start`, `end`) are the model's zero-based, half-open coordinates `[start, end)` (never silently converted to GenBank one-based display coordinates).

## ClusterBlast sidecar enrichment

- Parse text sidecars from canonical `clusterblast/`, `knownclusterblast/`, and `subclusterblast/` directories using natural region sorting (`_cN.txt`).
- Parse native antiSMASH JSON `antismash.modules.clusterblast` modules validating module schema `2` and result schema `5`.
- Precedence is per `(record_id, region_number, search_type)` key: text files are preferred when present for a key, while JSON fills remaining keys.
- Retain valid empty results as negative evidence and retain source path, SHA-256, source format, and supported schema versions as provenance.
- In strict mode, malformed or unattached sidecars raise `ClusterBlastParseError`. In lenient mode, errors emit a `clusterblast_parse_failed` diagnostic.
- Retain `misc_feature` entries only as raw GenBank evidence; do not infer that they are ClusterBlast results.

## Comparative review

- Compare two antiSMASH runs or files with `compare <left> <right>`.
- Supported record matching modes:
  - `record_id` (default): Exact record ID matching; requires unique IDs on each side.
  - `record_region`: Key by `(record_id, region_number)`; valid only when every record contains exactly one numbered region.
  - `single_record`: Explicit user-requested pairing of single-record inputs with differing IDs.
  - `coordinate_overlap`: Match records by reciprocal overlap of feature spans. Requires explicit `--assume-shared-coordinate-system` and a reciprocal overlap fraction in `(0, 1]` (default 0.80).
- The coordinate assumption is appropriate only when coordinate correspondence has been independently established, such as re-annotations of the same contigs. It does not establish sequence homology between arbitrary isolates.
- antiSMASH region GBKs commonly rebase extracted regions to coordinate zero. Such files do not automatically share their original chromosome coordinate system and must not use coordinate matching solely because they came from the same run.
- Comparison evaluates product gains/losses with multiset counts, new/resolved diagnostics, feature counts delta, and intergenic distance summaries (including circular topology wrap gaps).
- Comparative JSON exports `schema_name: "antismash-review-comparison"` with `schema_version: "0.1.0"`.
