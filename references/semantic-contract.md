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
- Preserve repeated `/specificity` values in order in both `Domain.specificity` and the raw
  qualifier map. The current model does not classify T2PKS/KR/Minowa semantics; downstream
  code must interpret the retained strings conservatively.
- Keep `aSModule` objects separate from domains. Resolve their domain IDs and diagnose missing or duplicated references.
- Allow `CDS_motif` records without a motif label or e-value, including RiPP prepeptide annotations.
- Retain raw `PFAM_domain` hits and expose a stable deduplicated view keyed by record, locus, nucleotide/protein coordinates, and version-normalized accession.

## Review boundaries

- `region.contig_edge=True` means the extracted BGC context reaches a record boundary. It does not prove that a biosynthetic core or CDS is truncated.
- Diagnose a proto-core at a record edge separately.
- Diagnose a CDS edge truncation only when its location is fuzzy/partial and reaches the edge.
- Do not emit a generic motif-confidence warning based on one e-value cutoff.
- Standalone `gene` features remain raw evidence. A `/pseudo` gene overlapping a region
  emits `pseudogene_in_cluster`; this does not prove a frameshift or functional loss.
- Feature types outside the adapter set are retained in `raw_features` and produce one
  aggregated `unrecognized_feature_type` NOTICE per record, except the structural `source`
  feature.

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
| `pseudogene_in_cluster` | WARNING | A standalone `/gene` with `/pseudo` overlaps one or more regions | Does not prove a frameshift, functional loss, or compound identity |
| `unrecognized_feature_type` | NOTICE | A non-structural feature type is retained only in `raw_features` | Does not mean the source annotation is invalid; it marks adapter coverage |
| `clusterblast_parse_failed` | WARNING (lenient) | A recognized ClusterBlast text or JSON sidecar could not be parsed or attached safely | Does not repair or silently discard the malformed sidecar; strict mode fails instead |

Region-level Pfam duplication does not create a diagnostic. Raw and deduplicated Pfam totals already expose the underlying representation.

## Discovery and exports

- Classify recognized files deterministically as native JSON, aggregate GenBank, or region GenBank.
- During inspection, consume aggregate GenBank files in preference to region files to avoid duplicate representations.
- `loading.py` owns the shared GenBank plus ClusterBlast enrichment path used by inspection and comparison; the compatibility loader in `cli.py` is only a wrapper.
- Native antiSMASH JSON alone remains unsupported as primary input; when a GenBank input is provided alongside JSON or sidecar text directories, sidecar enrichment parses and attaches ClusterBlast results.
- Markdown and the top-level JSON diagnostics include evidence-scoped review diagnostics.
- The JSON envelope reports `schema_name: "antismash-review"`, `schema_version: "0.3.0"`, and `parser_version`, and embeds `antismash_provenance` on each record.
- TSV (`--format tsv`) is a compact one-row-per-record summary.
- Entity-level TSV exports (`--format gene-tsv`, `--format domain-tsv`, and `--format clusterblast-tsv`) export one row per entity. All entity TSV coordinates (`start`, `end`) are the model's zero-based, half-open coordinates `[start, end)` (never silently converted to GenBank one-based display coordinates).

## Assembly-line evidence

- `predict_assembly_lines(record)` is a derived interpretation layer; it does not change the parsed `Record` JSON model.
- A local prediction groups modules only by explicit antiSMASH locus-tag membership or one explicit multi-CDS module. Separate CDS-local chains in one region remain separate outputs.
- Within one CDS, protein/domain coordinates are preferred for ordering. If they are unavailable, nucleotide coordinates are ordered with reverse-strand orientation. Cross-CDS proximity is never used as a Phase 1 ordering heuristic.
- `Module.monomer_pairings` have priority and are parsed as exact substrate-to-monomer evidence. Multiple calls, `X`, malformed calls, and non-proteinogenic tokens remain visible in the typed result and exports.
- When pairings are absent, only `specificity="substrate consensus: ..."` values on AMP-binding/A-like domains are interpreted as low-confidence fallback calls. Other specificity values remain raw notes; conflicting consensus values remain ambiguous.
- Thioesterase-like evidence is reported as release-domain evidence with release mode unknown. A terminal or final module flag does not prove hydrolysis or cyclization.
- Assembly-line mass values are separate modeled candidates, not fields on the parsed evidence objects. The dedicated `antismash-review-assemblyline` JSON schema is version `0.2.0` after the Phase 2 chemistry gate.
- Core-mass candidates are formula-derived only for fully resolved, high-confidence proteinogenic amino-acid-like NRPS chains using genuine free amino-acid formulas with peptide dehydration. Unknown `X`, non-proteinogenic calls, PKS modules, starter/acyl uncertainty, iterative modules, low-confidence specificity fallbacks, and tailoring chemistry make the full-core candidates null.
- When chemistry is complete enough, linear and head-to-tail cyclic candidate values are emitted separately with `topology_assumption="unknown"`; neither value is the measured final metabolite mass.

## Architecture assessment

- `assess_architecture(record)` returns scoped assessments rather than one record-wide score. T1PKS/trans-AT expectations are region-scoped; NRPS expectations are module-scoped.
- Canonical T1PKS requires parsed KS, AT, and ACP/PCP evidence. A trans-AT product label requires KS and ACP/PCP but explicitly exempts missing cis-AT evidence.
- A canonical NRPS starter module requires A and ACP/PCP; non-starter modules additionally require C. Final, iterative, and incomplete flags remain antiSMASH evidence and are not silently re-derived.
- Unsupported labels such as `NRPS-like`, T2PKS, T3PKS, prodigiosin, and other specialized classes return `not_applicable` until a class-specific rule is audited.
- Product labels beyond the canonical trio (T1PKS, transAT-PKS, NRPS) intentionally return `not_applicable`. This includes `NRPS-like`, `T2PKS`, `T3PKS`, `prodigiosin`, `lanthipeptide`, and all other specialized classes. Class-specific architecture scoring is not performed for these labels. Raw domain/module evidence remains available for qualified downstream interpretation by researchers or AI-assisted workflows.
- Architecture scores measure expected parsed-domain slot coverage, not probability of pathway completeness, activity, or metabolite production. The legacy `missing_nrps_pks_architecture` diagnostic remains available alongside the more specific missing-slot warning.

## Genome-browser tracks

- `review_findings(record)` preserves the existing diagnostics while attaching a structured location/entity when the rule has one. Exporters never parse diagnostic message text to rediscover coordinates.
- GFF3 converts internal `[start, end)` coordinates to one-based inclusive coordinates by adding one only to the start. BED retains internal zero-based half-open coordinates.
- GFF3 and BED export regions, collection layers, CDSs, domains, modules, and localized review findings. Compound locations produce one row per part with deterministic part IDs; cross-origin status is retained.
- Aggregate records use their record ID as `seqid`. A rebased region GenBank uses a source-qualified `seqid` so repeated `contig_1` labels are not merged. Duplicate preferred IDs receive deterministic suffixes.
- GFF3 emits `.` for score and phase; the parser has no validated CDS phase contract. Attributes are percent-encoded, and repeated renders are byte-identical.

## ClusterBlast sidecar enrichment

- Parse text sidecars from canonical `clusterblast/`, `knownclusterblast/`, and `subclusterblast/` directories using natural region sorting (`_cN.txt`).
- The current contract reads recognized `.txt` sidecars directly inside those canonical
  directories; nested HTML/detail assets are ignored.
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
  - `coordinate_overlap`: Match records by reciprocal overlap of feature spans. Requires explicit `--assume-shared-coordinate-system`; the optional `--min-reciprocal-overlap` override accepts a reciprocal overlap fraction in `(0, 1]` (default 0.80).
- The coordinate assumption is appropriate only when coordinate correspondence has been independently established, such as re-annotations of the same contigs. It does not establish sequence homology between arbitrary isolates.
- antiSMASH region GBKs commonly rebase extracted regions to coordinate zero. Such files do not automatically share their original chromosome coordinate system and must not use coordinate matching solely because they came from the same run.
- Comparison evaluates product gains/losses with multiset counts, new/resolved diagnostics, feature counts delta, and intergenic distance summaries (including circular topology wrap gaps).
- Comparative JSON exports `schema_name: "antismash-review-comparison"` with `schema_version: "0.2.0"`.

## Provenance

- Parsed records retain all antiSMASH structured-comment keys and repeated values in a typed internal provenance object. Biopython's structured-comment mapping is supplemented with the raw GenBank comment block when repeated keys would otherwise be lost. The existing record JSON envelope omits that object until a deliberate record-schema migration.
- The dedicated provenance manifest reports source paths, SHA-256 hashes, record IDs, review-tool version, normalized antiSMASH version/run date fields when present, database-like fields, and raw unknown keys.
- Missing metadata remain unknown. A comparison provenance delta uses `True`/`False` only when both values are present; otherwise it uses `None` and never calls two absent fields unchanged.
- Comparison schema version is `0.2.0` because matched record comparisons now carry an explicit provenance delta.

## Cohort matrices and clustering

- `cohort ROOT` treats each immediate child directory as one member and sorts member names deterministically. `cohort --manifest samples.tsv` accepts `sample<TAB>path` rows, preserves their listed order, and requires exactly one root or manifest.
- Each member uses the shared discovery/loading path. Aggregate GenBank input takes precedence over region GenBank input for counting, while all discovered source and sidecar paths remain protected from output overwrite.
- Product matrix counts are region-level product labels across all records in a member. Domain matrix counts are adapted `aSDomain` features only; `PFAM_domain` hits are not mixed into the domain matrix.
- Matrix keys use Unicode NFKC normalization, stripping, and case folding. JSON retains normalized keys plus canonical display and raw labels. Binary values mean member-level presence; count values are integer member totals.
- Invalid members fail by default with their name and path. `--skip-invalid-members` is an explicit mode and lists skipped members/errors in JSON; no invalid strain disappears silently.
- Cohort JSON uses `schema_name: "antismash-review-cohort"` and `schema_version: "0.1.0"`, including source hashes, provenance summaries, matrices, and null clustering fields when clustering is not requested.
- `--cluster-by domain-jaccard` computes binary-domain Jaccard distances, deterministic average-linkage clustering with lexical tie-breaking, and a valid Newick tree. Clustering is optional and adds no runtime dependency.
