# antiSMASH-review sixth extension plan

**Repository:** `WhyAdr/antiSMASH-review`
**Plan target:** repository state inspected at `main` commit `f112dd90541bd018e9496470a39dcdf76471634e` (`feat: implement synthetic test fixtures, feature diagnostics, and sidecar validation`, 2026-08-09)
**Current package version:** `0.1.0`
**Current record JSON schema:** `antismash-review` `0.2.0`
**Current comparison JSON schema:** `antismash-review-comparison` `0.2.0`

## Executive summary

This sixth pass should extend `antiSMASH-review` from a careful parser/reviewer into a modest **BGC interpretation and interoperability layer** without blurring the line between antiSMASH evidence and downstream inference.

The proposed work has five feature areas, split into seven implementation phases:

1. **NRPS/PKS assembly-line interpretation** — reconstruct antiSMASH-predicted module/monomer chains, preserve ambiguity, and identify likely release domains. Conservative core-scaffold mass is a separate Phase 2 gate.
2. **Domain-composition-aware architecture assessment** — replace the current binary “NRPS/PKS product but no NRPS/PKS domains” check with product- and module-aware minimum-domain expectations, including explicit exceptions for trans-AT PKSs, starter/final modules, iterative modules, and incomplete/edge-truncated loci.
3. **N-way cohort analysis** — generalize reusable feature fingerprints beyond the current two-input comparison and add a `cohort` workflow for product-class and domain-content matrices, optional deterministic domain-Jaccard clustering, and machine-readable cohort output.
4. **GFF3/BED interoperability** — export antiSMASH collections, CDSs, domains, modules, and localized review findings as genome-browser tracks while preserving the package's existing zero-based internal coordinate contract.
5. **Run provenance manifest** — capture antiSMASH run metadata from the structured comment, preserve unknown fields, expose source hashes and parser versions, and make provenance differences visible during run-to-run comparison.

The main architectural recommendation is **not** to put all of this logic into `review.py`, `compare.py`, or `cli.py`. The parser models are already strong enough to support these additions. Instead, add small analysis modules with typed result objects, then let review rules, exporters, `compare`, and the new `cohort` command consume the same derived primitives.

The original feature-area order is:

```text
Area A  NRPS/PKS assembly-line prediction
Area B  domain-aware architecture assessment
Area C  GFF3/BED interoperability
Area D  provenance manifest
Area E  cohort matrices and optional clustering
```

`cohort` is intentionally last because it is the largest surface-area change. The first four phases are individually useful for day-to-day BGC review and bench follow-up.

## Canonical implementation split

This is an umbrella roadmap, not one implementation patch or one release. Each phase below
must have its own commit or small commit group, tests, and acceptance gate.

```text
Phase 0  baseline, fixture policy, shared loader, fingerprints
Phase 1  evidence-only NRPS/PKS assembly-line reconstruction
Phase 2  optional chemistry registry and conservative core-mass candidates
Phase 3  architecture assessment and diagnostics
Phase 4  localized findings plus GFF3/BED export
Phase 5  provenance survey, manifest export, and comparison deltas
Phase 6  cohort product/domain matrices
Phase 7  optional deterministic domain-Jaccard clustering
```

The required first deliverable is Phase 0 plus Phase 1. Phase 2 is a separate scientific
decision gate. If the chemistry registry cannot be validated from real antiSMASH calls and
independently calculated fixtures, stop after Phase 1 and ship evidence-only assembly-line
output. Clustering must not block the cohort matrix release.

The current execution has passed the Phase 0 through Phase 7 gates. The implementation keeps
provenance, cohort matrices, and optional clustering behind their own schemas and CLI contracts.

The five detailed topic sections below are design specifications mapped onto these gates:

| Detailed section | Canonical phase | Release boundary |
|---|---:|---|
| Shared loader and fingerprints | 0 | No serialized schema change |
| Assembly-line reconstruction | 1 | Derived assembly-line evidence schema |
| Chemistry and core mass | 2 | Separate opt-in fields/export; null for unresolved chemistry; assembly-line schema 0.2.0 |
| Architecture assessment | 3 | New diagnostics only after fixture-backed rule review |
| Localized findings and GFF3/BED | 4 | New exporters; existing review API remains compatible |
| Provenance | 5 | Manifest first; `Record` schema bump only if embedding provenance |
| Cohort matrices | 6 | Product/domain matrices without clustering |
| Domain-Jaccard clustering | 7 | Optional follow-on; not required for cohort output |

---

# 1. Current-state findings that shape this plan

## 1.1 The parser already captures most of the information needed

The current typed model contains:

- `Domain.domain_id`
- `Domain.name`
- `Domain.subtypes`
- `Domain.tool`
- `Domain.locus_tag`
- `Domain.protein_start` / `protein_end`
- `Domain.specificity`
- `Module.domain_ids`
- `Module.locus_tags`
- `Module.module_type`
- `Module.complete`
- `Module.starter`
- `Module.final`
- `Module.iterative`
- `Module.monomer_pairings`
- `Module.multi_cds`
- `Module.missing_domain_ids`

This means the sixth pass does **not** need a second parser for NRPS/PKS logic. The missing layer is interpretation of already-parsed antiSMASH evidence.

The upstream antiSMASH `Module` implementation is especially important here. antiSMASH serializes substrate/monomer calls explicitly as:

```text
substrate -> monomer
```

and permits a module to span more than one CDS, specifically noting split-CDS trans-AT PKS cases as a legitimate use case. Therefore:

- `Module.monomer_pairings` should be the primary source for monomer calls;
- a multi-CDS module must **not** be automatically flagged as malformed;
- the reviewer should preserve antiSMASH's own `complete`, `starter`, `final`, and `iterative` flags rather than re-deriving them from scratch.

## 1.2 The current architecture diagnostic is intentionally coarse

`review.py` currently emits `missing_nrps_pks_architecture` only when a region product contains `NRPS` or `PKS` but the record has no domains whose tool is exactly `nrps_pks_domains`.

That is a useful sanity check, but it cannot distinguish:

- a canonical cis-AT T1PKS lacking AT evidence;
- a legitimate trans-AT PKS that should not contain a cis-AT in every module;
- an NRPS starter module that may lack a condensation domain;
- an explicitly incomplete antiSMASH module;
- a region-edge truncation;
- an NRPS/PKS-like product label for which a canonical module expectation would be overconfident.

The next step should therefore be a typed **architecture assessment** rather than more keyword branches directly in `review.py`.

## 1.3 Pairwise comparison is a good base, but not a cohort engine

`compare.py` currently compares exactly two record sets and supports:

- `record_id`
- `record_region`
- `single_record`
- guarded `coordinate_overlap`

It also already contains useful reusable ideas such as product multisets, diagnostic fingerprints, feature counts, and intergenic summaries.

For N-way analysis, do not repeatedly invoke pairwise comparison for every strain pair. Instead, extract shared feature/fingerprint functions into a small reusable module and build cohort matrices directly from each member's records.

## 1.4 Entity TSVs establish a useful exporter pattern

`exporters/entity_tables.py` already provides deterministic row-oriented exports and explicitly preserves the internal zero-based, half-open coordinate model. GFF3/BED should follow the same pattern rather than embedding formatting logic into the CLI.

Critical coordinate rule:

```text
internal model / BED:  zero-based, half-open  [start, end)
GFF3:                  one-based, inclusive   start+1 .. end
```

That conversion must be covered by explicit off-by-one tests.

## 1.5 Adding fields to `Record` changes the JSON schema automatically

`exporters/json_export.py` currently serializes records with `dataclasses.asdict(record)`.

Therefore, adding a new field such as `Record.provenance` is **not** an invisible internal change. It changes the public JSON envelope and must trigger a schema version bump under the repository's existing compatibility policy.

Recommendation: group all intended `Record` model additions for this pass and bump the record schema **once**, from `0.2.0` to `0.3.0`, rather than incrementing it in several small commits.

## 1.6 Baseline after the fifth pass

The fifth pass is implemented in `f112dd9`, not `af23c85`. The current public suite reports
57 passing tests; Ruff and strict mypy pass. The synthetic GenBank builder also emits
Biopython warnings for antiSMASH qualifier keys longer than the GenBank standard's nominal
20-character width. Treat warning cleanup or an explicit warning policy as a Phase 0 task;
do not let the sixth pass multiply warning noise.

The local `PA-LS101--GCA-030388725.1-antismash/` bundle is a useful private integration input
but is currently untracked and is not covered by the existing exact-directory ignore rule.
Keep it outside the distributable tree or add it to local Git excludes. Never make it a
required committed fixture. A read-only audit of its aggregate GenBank found two records; the
primary `contig_1` record contains 11 regions, 10 modules, a real multi-CDS module, reverse-
strand modules, incomplete/unknown module states, and monomer calls including `D-Orn`, `Ser`,
`His`, `X`, `ccmal`, and `NH2`. These observations should drive optional integration tests and
rule-registry surveys, not be copied into public synthetic fixtures.

## 1.7 Execution status after Phase 7

The Phase 0-7 implementation in this checkout reports 129 passing tests, clean Ruff and strict
mypy checks, deterministic cohort matrix and clustering fixtures, and a no-isolation wheel
validation gate. A read-only PA-LS101 validation also produced one provenance input, 13 product
columns in cohort JSON, and 12 assembly-line TSV rows. The private bundle remains ignored and is
not required for the test suite.

---

# 2. Phase 0 - baseline, shared loading, fingerprints, and contracts

Before adding analytical features, make two small structural changes. This phase must not add
assembly-line inference, architecture warnings, provenance fields, GFF/BED formats, or cohort
behavior.

### Phase 0 entry and exit gate

Entry: the fifth-pass baseline is green and private integration directories are ignored or
outside the checkout.

Exit:

- `loading.py` owns one discover -> parse -> sidecar-enrichment path;
- `compare` produces byte-equivalent output for the existing public comparison fixtures;
- fingerprint helpers have explicit raw-versus-normalized semantics;
- `DiagnosticFingerprint` has one owner, avoiding a `fingerprints.py` <-> `compare.py`
  import cycle;
- no existing JSON, TSV, matching-mode, or strict/lenient behavior changes; and
- Ruff, mypy, tests, coverage, and wheel/import checks pass.

## 2.1 Move record loading out of `cli.py`

Current public-ish helper:

```python
load_review_records(manifest, *, lenient) -> tuple[list[Record], set[Path]]
```

lives in `cli.py`. `compare` already depends on it indirectly through CLI orchestration, and `cohort` will need exactly the same discovery + GenBank + ClusterBlast enrichment behavior.

Move it to:

```text
antismash_review/loading.py
```

Suggested API:

```python
@dataclass(slots=True)
class LoadedReviewInput:
    root: Path
    records: list[Record]
    input_paths: set[Path]


def load_review_input(
    manifest: InputManifest,
    *,
    lenient: bool = False,
) -> LoadedReviewInput:
    ...
```

Keep a compatibility import or thin wrapper in `cli.py` during this release if tests or downstream code currently import `load_review_records`.

This gives `inspect`, `compare`, and `cohort` one canonical loading path.

## 2.2 Add reusable record fingerprints

Create:

```text
antismash_review/fingerprints.py
```

with dependency-free helpers such as:

```python
def product_counter(records: Sequence[Record]) -> Counter[str]: ...
def domain_counter(records: Sequence[Record]) -> Counter[str]: ...
def domain_presence(records: Sequence[Record]) -> frozenset[str]: ...
def diagnostic_counter(records: Sequence[Record]) -> Counter[DiagnosticFingerprint]: ...
```

Move `DiagnosticFingerprint` to `fingerprints.py` (or to `models.py`) before adding
`diagnostic_counter`; `fingerprints.py` must not import it back from `compare.py`. `compare.py`
should consume these helpers instead of maintaining its own product extraction logic.
`cohort.py` can then use the same semantics.

### Normalization rule

Be conservative. For product names and domain names:

- strip surrounding whitespace;
- preserve the original antiSMASH token for display;
- use a documented normalized key only for matrix identity, e.g. Unicode-normalized + case-folded;
- do not silently map biologically distinct product classes into one bucket.

If synonyms are introduced later, make the synonym registry explicit and export both raw and
normalized values. Preserve the existing pairwise comparison's raw product-delta behavior in
Phase 0; use normalized keys for new matrices only until a deliberate comparison-schema change
is approved.

## 2.3 Schema policy for this pass

Recommended schema changes:

```text
RECORD_SCHEMA_VERSION      0.2.0 -> 0.3.0   # only if Record gains provenance
COMPARISON_SCHEMA_VERSION  0.1.0 -> 0.2.0   # if provenance deltas are embedded
COHORT_SCHEMA_NAME         antismash-review-cohort
COHORT_SCHEMA_VERSION      0.1.0
```

Assembly-line predictions can initially use a separate derived-result schema rather than being inserted into every `Record` JSON object. That reduces coupling between parser evidence and analytical inference.

Suggested new schema:

```text
ASSEMBLYLINE_SCHEMA_NAME     antismash-review-assemblyline
ASSEMBLYLINE_SCHEMA_VERSION  0.2.0  # after Phase 2 core-mass candidates
```

---

# 3. Phase 1 - evidence-only NRPS/PKS assembly-line reconstruction

## 3.1 Goal

Turn the existing module/domain evidence into a human- and machine-readable hypothesis such as:

```text
Assembly line 1
  M1  Ser
  M2  Leu
  M3  Thr
  release: thioesterase-like domain present; release mode unresolved

Predicted core chain: Ser–Leu–Thr
Core neutral monoisotopic mass: 319.1743 Da
Mass coverage: 3/3 monomers resolved
```

The central rule is:

> This predicts an antiSMASH-derived **core scaffold hypothesis**, not the final metabolite identity or measured molecular mass.

For lipopeptides, glycosylated products, heavily tailored metabolites, hybrid NRPS-PKS systems, unusual extender units, or unresolved monomers, the output must say exactly which chemistry is missing.

Phase 1 is the evidence-only assembly-line layer. Phase 2 adds a separate mass object only
when its chemistry gate passes; unresolved or unsupported chemistry continues to emit null
full-core candidates and explicit caveats.

## 3.2 New module

Create:

```text
antismash_review/assemblyline.py
```

Keep this separate from `review.py`. Parsing should remain evidence capture; assembly-line interpretation is derived analysis.

Suggested dataclasses:

```python
@dataclass(slots=True, frozen=True)
class MonomerCall:
    substrate: str | None
    monomer: str | None
    display: str
    source: Literal["module_pairing", "domain_specificity", "unknown"]
    confidence: Literal["high", "medium", "low", "unresolved"]
    notes: tuple[str, ...] = ()


@dataclass(slots=True, frozen=True)
class ModulePrediction:
    index: int
    module_type: str | None
    locus_tags: tuple[str, ...]
    domain_ids: tuple[str, ...]
    domain_names: tuple[str, ...]
    complete: bool
    starter: bool
    final: bool
    iterative: bool
    multi_cds: bool
    monomer_calls: tuple[MonomerCall, ...]
    release_domains: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(slots=True, frozen=True)
class MassEstimate:
    linear_core_mass_da: float | None
    head_to_tail_cyclic_candidate_mass_da: float | None
    resolved_monomers: int
    total_monomers: int
    coverage_fraction: float
    unresolved_monomers: tuple[str, ...]
    topology_assumption: Literal["linear", "cyclic", "unknown"]
    chemistry_scope: str


@dataclass(slots=True, frozen=True)
class AssemblyLinePrediction:
    record_id: str
    region_number: int | None
    strand: int | None
    modules: tuple[ModulePrediction, ...]
    chain: tuple[MonomerCall, ...]
    mass: MassEstimate | None = None
    ordering_basis: str
    ordering_confidence: Literal["high", "medium", "low"]
    caveats: tuple[str, ...]
```

Public function:

```python
def predict_assembly_lines(record: Record) -> list[AssemblyLinePrediction]:
    ...
```

### Phase 1 entry and exit gate

Entry: Phase 0 has frozen loader/fingerprint behavior and a fixture survey has established the
actual module/domain qualifier forms used by the supported antiSMASH output.

Exit:

- `Module.monomer_pairings` are parsed without discarding raw strings;
- `X`, missing pairings, conflicting calls, non-proteinogenic calls, and multi-CDS modules
  remain explicit rather than being silently normalized;
- module grouping is deterministic and does not merge separate CDS-local chains merely because
  they share a region;
- reverse-strand order is tested using protein coordinates or strand-aware nucleotide order;
- no mass calculation is required for ordering; Phase 2 mass candidates remain a separate
  derived object and are null when chemistry is unsupported; and
- a derived assembly-line export can be rendered without changing `Record` JSON.

## 3.3 Resolve each module's domains once

`genbank._resolve_modules()` currently checks that domain IDs exist, but does not expose the resolved objects.

Do not duplicate ad hoc lookup code in several features. Add a small helper:

```python
def domains_for_module(record: Record, module: Module) -> list[Domain]:
    ...
```

It should:

- retain `module.domain_ids` order;
- return only uniquely resolved domains;
- rely on existing `module_domain_missing` / `domain_id_duplicated` diagnostics for malformed references;
- never guess which duplicated domain is canonical.

A later optimization can cache the ID index, but records are small enough that correctness matters more than micro-optimization.

## 3.4 Monomer-call priority

Use a strict evidence hierarchy.

### Priority 1 — `Module.monomer_pairings`

This is already antiSMASH's assembled substrate-to-monomer result. Parse exact strings of the form:

```text
substrate -> monomer
```

The parser currently preserves these strings verbatim, so the analysis layer should validate rather than assume they are always well formed.

If more than one pairing occurs in a module, retain all calls and mark the module as ambiguous/multi-call instead of arbitrarily choosing one.

### Priority 2 — recognized A-domain specificity consensus

Only use `Domain.specificity` as a fallback when the module has no monomer pairing.

Do **not** simply take the first specificity string. antiSMASH may preserve multiple predictor outputs. Add a dedicated parser for recognized consensus forms after inspecting fixtures from the supported antiSMASH versions.

Until a specificity form is explicitly supported:

```text
source = domain_specificity
confidence = low/unresolved
```

and preserve the raw strings in notes.

### Priority 3 — unresolved

If neither evidence source yields a defensible monomer:

```text
M4 = ?
```

This is preferable to converting uncertainty into a false sequence.

## 3.5 Assembly-line ordering

This is the scientifically delicate part.

### Do not assume that all modules in one BGC form one linear chain

A region can contain:

- multiple NRPS/PKS proteins;
- parallel or branching assembly lines;
- separate biosynthetic systems inside a large hybrid region;
- modules on opposing strands;
- trans-acting enzymes;
- split-CDS modules;
- trans-AT PKSs whose genomic organization is not equivalent to a simple cis-AT peptide-like assembly line.

### Recommended ordering tiers

Phase 1 executes Tier A and Tier B only. Tier C is retained below as a deferred design note,
not as an implementation requirement.

#### Tier A — within one CDS: high confidence

For modules that share one locus tag:

- sort by protein/domain order when protein coordinates are available;
- otherwise sort by nucleotide coordinates respecting strand.

For a negative-strand CDS, biological order is reverse genomic coordinate order.

#### Tier B — one antiSMASH multi-CDS module: high confidence for module membership

Treat the module as one module because antiSMASH explicitly encoded it that way. Preserve the ordered `locus_tags`/domain evidence; do not split it merely because `multi_cds=True`.

#### Tier C - across multiple CDSs: deferred in Phase 1

Phase 1 must not combine separate CDS-local module chains using a proximity heuristic. The
current model does not contain an authoritative assembly-line graph, and a shared region can
contain parallel, hybrid, or unrelated systems. Emit one local chain per explicit antiSMASH
module group and preserve separate CDS-local chains as separate outputs.

If a later phase introduces cross-CDS candidates, it must expose the candidate member modules,
ordering basis, competing alternatives, and `ordering_confidence="medium"` or `"low"`. It must
never replace multiple plausible chains with one canonical sequence.

The following criteria are deferred design inputs for that later phase:

- same strand;
- adjacent or near-adjacent core biosynthetic CDSs;
- same protocluster/region;
- compatible module types;
- no competing alternative ordering.

Label this:

```text
ordering_confidence = medium
ordering_basis = genomic-coorientation
```

If ordering is ambiguous, emit multiple local chains rather than one false canonical product sequence.

A future extension can use antiSMASH JSON module ordering or graph-level assembly-line information if a stable upstream representation exists.

## 3.6 Release-domain handling

Claude's illustrative `TE(cyclize)` is useful as a display concept but too specific as a default biological claim.

A thioesterase-like terminal domain can support chain release by hydrolysis or macrocyclization depending on the system. Unless antiSMASH provides an explicit release-mode call, render:

```text
TE(release; mode unknown)
```

not:

```text
TE(cyclize)
```

Likewise, a `final_module` flag is evidence that antiSMASH considers the module terminal; it is not itself proof of a particular chemical termination reaction.

Add a normalized release-domain alias registry populated from observed antiSMASH `aSDomain` names and tested fixtures.

## 3.7 Conservative core-mass calculation

This is a separate Phase 2 deliverable and must not be implemented in the Phase 1
assembly-line commit.

### Phase 2 entry and exit gate

Entry: Phase 1 has only evidence-level monomer calls and a fixture table of exact observed
monomer strings. No chemistry alias is accepted from an unverified example.

Exit:

- free-monomer versus residue-mass semantics are explicit;
- linear and cyclic candidates are separate fields;
- unresolved `X`, `NH2`, `ccmal`, unknown tails, hybrid PKS chemistry, and tailoring return
  null full-core candidates with machine-readable caveats;
- independently calculated fixture masses pass within a documented tolerance; and
- mass output is opt-in and does not append columns to existing stable TSV contracts.

### Scope for v1

Mass estimation should initially support **fully resolved amino-acid-like NRPS cores**.

Do not calculate a final “natural product mass” from:

- unresolved monomers;
- unknown fatty-acyl starters;
- PKS extender units whose reduction state is not explicitly modeled;
- glycosylation;
- halogenation;
- oxidation/reduction not represented in the monomer call;
- macrocyclization unless topology is explicitly known;
- other tailoring enzymes outside the module call.

For hybrid NRPS-PKS chains, emit the chain but return a mass only if every incorporated monomer has a curated mass model.

### Mass registry

Create:

```text
antismash_review/chemistry.py
```

with a small curated table rather than pulling in a cheminformatics dependency.

Start with:

- the 20 proteinogenic amino acids;
- exact aliases observed in antiSMASH monomer outputs;
- stereochemical aliases where mass is unchanged, e.g. L/D forms;
- additional non-proteinogenic monomers only when a fixture and a literature/chemical formula source are documented in code comments/tests.

Do not attempt a free-form parser that interprets every `Me`, `OH`, `D-`, or unusual abbreviation in the first implementation. Exact known aliases are safer.

### Formula

For a linear peptide with `n` free amino-acid monomers, where each registry value is a
free-monomer monoisotopic mass:

```text
M_linear = sum(M_free_monomer) - (n - 1) * M_H2O
```

For a head-to-tail cyclic peptide formed by one additional condensation:

```text
M_cyclic = M_linear - M_H2O
```

However, because release topology is often unresolved, the default should be:

```text
topology_assumption = unknown
```

and the primary reported value should be the linear core mass only when that assumption is explicitly stated.

A safer alternative is to report both candidate values when every monomer is known:

```text
linear_core_mass_da
head_to_tail_cyclic_candidate_mass_da
```

with neither labeled as the observed metabolite mass. Do not overload one generic mass field: a
linear candidate, a cyclic candidate, and an observed metabolite mass are different claims.

### Required uncertainty metadata

Every mass result must include:

- resolved monomer count;
- total monomer count;
- coverage fraction;
- unresolved monomer names;
- whether a fatty-acyl starter is unresolved;
- whether a PKS/hybrid segment prevents complete calculation;
- topology assumption;
- explicit chemistry scope.

If coverage is less than 100%, default behavior should be:

```text
linear_core_mass_da = null
head_to_tail_cyclic_candidate_mass_da = null
```

Optionally expose a separately named `resolved_partial_mass_da`; never present a partial value as the full scaffold mass.

## 3.8 Export surfaces

Add:

```text
antismash_review/exporters/assemblyline_markdown.py
antismash_review/exporters/assemblyline_json.py
antismash_review/exporters/assemblyline_table.py
```

Suggested CLI formats:

```bash
antismash-review inspect result/ --format assemblyline-tsv
antismash-review inspect result/ --format assemblyline-json
```

The normal Markdown `inspect` report may include a concise “Predicted assembly lines” section, but keep the detailed machine representation in dedicated exports.

Example TSV columns:

```text
record_id
region_number
assembly_line
module_index
locus_tags
module_type
complete
starter
final
iterative
multi_cds
domain_names
substrate
monomer
call_source
call_confidence
release_domains
linear_core_mass_da
mass_coverage
unresolved_components
```

## 3.9 Acceptance tests

Add synthetic fixtures/tests for:

Phase 1 must pass items 1-8 and must not assert a mass value. Phase 2 owns items 9-12 and
requires an independent chemistry fixture before those tests are enabled.

1. three-module forward-strand NRPS with complete monomer pairings;
2. reverse-strand NRPS ordering;
3. multi-CDS module that remains one module;
4. two distinct NRPS chains in one region that must not be collapsed;
5. missing monomer pairing with recognized specificity fallback;
6. conflicting specificity outputs retained as ambiguity;
7. terminal TE-like domain rendered as release-mode unknown;
8. starter/final/iterative flags preserved;
9. unknown monomer causes full mass to be `null`;
10. exact mass for a fully resolved canonical peptide matches an independently calculated fixture value within a documented tolerance;
11. lipopeptide-like starter produces an explicit unresolved-acyl-tail caveat;
12. hybrid PKS segment prevents peptide-only mass from being mislabeled as full scaffold mass.

---

# 4. Phase 3 - domain-composition-aware architecture assessment

## 4.1 Goal

Do not remove or silently change the existing coarse warning in this phase. Add a typed
assessment alongside it, then add more specific diagnostics only after the rule matrix and
false-positive fixtures are accepted.

### Phase 3 entry and exit gate

Entry: Phase 1 has deterministic local assembly-line output, and the private/public fixture
survey has established the actual `aSDomain` names and product labels to classify.

Exit:

- `DomainSlot` and every alias registry entry have source-verified examples;
- assessment scope is explicit (region, protocluster, or module) and never compares unrelated
  domains across the whole record;
- trans-AT, starter, final, iterative, incomplete, and edge cases have dedicated fixtures;
- the legacy `missing_nrps_pks_architecture` diagnostic remains compatible for one cycle;
- detailed assessments have a documented Python/JSON export surface; and
- a new warning is not enabled merely because a slot is absent from a synthetic record.

Replace the coarse review-only decision with a typed assessment, but retain the existing
diagnostic during the compatibility cycle. The analysis should be able to say:

```text
T1PKS architecture: partial
expected core slots: KS, AT, ACP/PP-binding
observed: KS, ACP/PP-binding
missing: AT
interpretation: cis-AT core expectation not met; possible fragmented annotation,
                atypical architecture, or incorrect product classification
```

while correctly handling:

```text
transAT-PKS: absence of a cis-AT domain is not itself an error
```

## 4.2 Separate assessment from diagnostics

Create:

```text
antismash_review/architecture.py
```

Suggested models:

```python
@dataclass(slots=True, frozen=True)
class ArchitectureExpectation:
    name: str
    product_keys: tuple[str, ...]
    required_slots: tuple[DomainSlot, ...]
    optional_slots: tuple[DomainSlot, ...] = ()
    notes: tuple[str, ...] = ()


@dataclass(slots=True, frozen=True)
class DomainSlot:
    name: str
    aliases: tuple[str, ...]
    required: bool = True
    evidence_scope: str = "module"


@dataclass(slots=True, frozen=True)
class ArchitectureAssessment:
    product: str
    status: Literal["complete", "partial", "ambiguous", "not_applicable"]
    score: float | None
    expected_slots: tuple[str, ...]
    observed_slots: tuple[str, ...]
    missing_slots: tuple[str, ...]
    exemptions: tuple[str, ...]
    evidence_domains: tuple[str, ...]
    caveats: tuple[str, ...]
```

`review.py` should turn assessments into a small set of diagnostics; the detailed assessment belongs in the analysis/export layer.

Define the public analysis boundary explicitly:

```python
def assess_architecture(record: Record) -> list[ArchitectureAssessment]:
    ...
```

Assess each region/module in its own scope. Do not calculate one record-wide slot score from
domains belonging to unrelated regions. Add an `architecture-json` or equivalent derived
export before exposing the assessment through the CLI; a diagnostic-only implementation would
discard the very detail needed to audit a warning.

## 4.3 Domain-family normalization

Create one explicit alias registry, e.g.:

```python
DOMAIN_FAMILIES = {
    "KS": {...},
    "AT": {...},
    "ACP_PCP": {...},
    "A": {...},
    "C": {...},
    "TE": {...},
}
```

Populate aliases from:

- current antiSMASH fixtures;
- upstream antiSMASH NRPS/PKS domain names;
- the current parser's `Domain.name` / `subtypes` values.

Do not identify core NRPS/PKS architecture from arbitrary Pfam descriptions if the corresponding `aSDomain` call is available. Pfam can be supporting evidence later, but mixing detection layers immediately would make rule semantics harder to audit.

## 4.4 Seed only conservative product rules

Start with product labels whose expected architecture is strong enough to audit.

### Canonical T1PKS

At the module/core level, expect evidence for:

```text
KS + AT + ACP/PP-binding
```

but use antiSMASH's own module boundaries where possible.

Important: upstream antiSMASH detection of `T1PKS` is based on a CDS containing a PKS AT and a ketosynthase-family hit. The review layer is performing a **stricter completeness audit**, not recreating the antiSMASH detection rule.

### transAT-PKS

Do **not** require cis-AT in every module.

Minimum module expectation should focus on the trans-AT architecture actually represented by antiSMASH, such as KS + carrier-protein evidence, with trans-acting AT support assessed separately at the region level.

This class needs an explicit exemption:

```text
missing cis-AT is expected/allowed for transAT-PKS
```

### NRPS

For a canonical elongation module, expect:

```text
A + PCP/PP-binding
```

with condensation-domain expectation conditional on module position/type.

Examples:

- starter module: lack of an upstream C domain may be normal;
- internal elongation module: C + A + PCP is the canonical expectation;
- final module: termination domain may be present but is not universally mandatory;
- iterative module: repeated use means module count is not equal to product residue count.

### Classes to defer initially

Do not force canonical module rules onto:

- `NRPS-like`
- `PKS-like`
- T2PKS
- T3PKS
- prodigiosin/nontraditional PKS
- PUFA synthases
- hglE-KS
- other specialized PKS classes

until class-specific expectations are defined.

For these, return:

```text
status = not_applicable
```

rather than a false incompleteness score.

## 4.5 Scoring semantics

A numeric score is useful for sorting but dangerous if interpreted as a probability of functionality.

Use a transparent slot-coverage score:

```text
score = satisfied_required_slots / total_required_slots
```

with explicit exemptions removed from the denominator.

Example:

```text
T1PKS expected: KS, AT, ACP
observed: KS, ACP
score = 2/3 = 0.667
status = partial
```

Document clearly:

> The architecture score measures expected parsed-domain coverage. It is not a probability that the BGC is complete, active, or capable of producing a metabolite.

## 4.6 Region/module context

Architecture must be assessed in the correct genomic context.

Add shared location helpers rather than manually repeating overlap logic from `genbank.py`:

```text
antismash_review/locations.py
```

Possible functions:

```python
def overlaps(left: Location, right: Location) -> bool: ...
def containing_regions(record: Record, location: Location) -> list[CollectionFeature]: ...
def domains_in_collection(record: Record, collection: CollectionFeature) -> list[Domain]: ...
def modules_in_collection(record: Record, collection: CollectionFeature) -> list[Module]: ...
```

Refactor `_overlaps()` in `genbank.py` to use the shared helper.

This becomes useful for GFF3 diagnostics and cohort region summaries too.

## 4.7 Interaction with antiSMASH completeness flags and contig edges

If antiSMASH already marks a module `complete=False`, expose that fact as evidence.

Do not override it with a home-grown binary label.

If a missing core domain occurs in a region with `contig_edge=True` or a partial CDS at the relevant edge, add a caveat such as:

```text
missing core-domain evidence may reflect record/assembly truncation
```

Do not automatically downgrade every edge-associated warning to a notice; a missing core domain remains important. The edge context changes interpretation, not the observation.

## 4.8 Review diagnostics

Retain the existing `missing_nrps_pks_architecture` for the extreme case of no parsed NRPS/PKS domains, at least for one compatibility cycle.

Add a more specific diagnostic, for example:

```text
architecture_core_domain_missing
```

Message template:

```text
T1PKS region 4 lacks expected core domain slot(s): AT.
Observed slots: KS, ACP/PP-binding. This is an architecture-consistency
warning, not evidence that the locus is nonfunctional.
```

Potential second diagnostic:

```text
architecture_module_incomplete
```

only if there is value beyond antiSMASH's own `Module.complete=False` flag. Avoid duplicating upstream information as noise.

## 4.9 Tests

At minimum:

- T1PKS KS-AT-ACP: complete;
- T1PKS KS-ACP: missing AT warning;
- transAT-PKS KS-ACP without cis AT: no false missing-AT warning;
- NRPS starter A-PCP: allowed starter architecture;
- internal NRPS C-A-PCP: complete;
- internal NRPS A-PCP: partial/missing C;
- explicitly incomplete module retains antiSMASH state;
- iterative module does not multiply residue count automatically;
- unsupported product class returns `not_applicable`;
- region-edge missing domain carries truncation caveat;
- no NRPS/PKS domains still triggers the existing coarse diagnostic.

---

# 5. Phase 4 - GFF3 and BED export

## 5.1 Goal

Make review evidence directly viewable alongside ONT alignments, coverage, variants, or other genome-browser tracks.

### Phase 4 entry and exit gate

Entry: `review_findings()` has structured locations for the findings selected for export, and
the input loader identifies whether coordinates came from an aggregate record or a rebased
region GenBank.

Exit:

- `seqid` is deterministic and documents whether it names a full record or an extracted region;
- rebased region records with repeated `record_id` values cannot be rendered as one apparent
  chromosome;
- GFF3 and BED coordinate conversions pass exact `[0, 100)` and compound/cross-origin tests;
- duplicate preferred IDs receive deterministic disambiguating suffixes;
- GFF3 attribute escaping, strand, phase policy, and BED column semantics are documented; and
- output is byte-identical across repeated renders.

Example:

```bash
antismash-review inspect result.gbk --format gff3 --output review.gff3
antismash-review inspect result.gbk --format bed --output review.bed
```

## 5.2 Exporter modules

Add:

```text
antismash_review/exporters/gff3.py
antismash_review/exporters/bed.py
```

Keep coordinate conversion private to the exporter.

## 5.3 GFF3 coordinate contract

The model remains zero-based, half-open.

Convert only at serialization:

```python
gff_start = part.start + 1
gff_end = part.end
```

Add an explicit code comment:

```python
# Internal/Biopython coordinates are zero-based half-open [start, end).
# GFF3 is one-based inclusive, so only start is incremented.
```

This comment is worth keeping because an apparently “symmetrical” end conversion would create an off-by-one bug.

## 5.4 BED coordinate contract

BED uses the model's native convention:

```text
chromStart = start
chromEnd   = end
```

No coordinate conversion.

## 5.5 Feature types to export

Minimum:

- `region`
- `cand_cluster`
- `protocluster`
- `proto_core`
- `CDS`
- `aSDomain`
- `aSModule`

Although Claude's initial suggestion mentioned `Gene`, `Domain`, and `CollectionFeature`, modules are especially useful once the assembly-line functionality exists and cost little to include.

## 5.6 Stable IDs

GFF3 needs deterministic identifiers even when antiSMASH features lack a locus tag/domain ID.

Create one ID helper:

```python
def stable_feature_id(
    record_id: str,
    feature_type: str,
    ordinal: int,
    preferred: str | None,
) -> str:
    ...
```

Examples:

```text
contig_1:CDS:EBIMEA_021660
contig_1:aSDomain:AMP-binding.3
contig_1:region:11
contig_1:aSModule:0007
```

Fallback IDs must be deterministic from record/feature order, not Python object identity.
Preferred IDs are not assumed unique: duplicate domain IDs, locus tags, or repeated record
labels must receive deterministic suffixes within the rendered document.

Define `seqid` before implementing either exporter. For aggregate GenBank records, the default
may be the record ID. For rebased region GenBanks, use a stable region-qualified identifier
derived from source filename plus record ID (or require an explicit `--seqid-mode`), and retain
the original record ID as an attribute. Never merge extracted region coordinates solely because
they share `contig_1`.

## 5.7 Attributes

Useful GFF3 attributes include:

### Collections

```text
ID
Name
product
number
contig_edge
creating_tool
```

### CDS

```text
ID
locus_tag
gene
product
gene_kind
region_numbers
```

### Domain

```text
ID
Name
domain_id
tool
locus_tag
subtypes
specificity
```

### Module

```text
ID
type
locus_tags
complete
starter
final
iterative
monomer_pairings
```

Use proper GFF3 percent-encoding for attribute values. Do not manually concatenate unescaped free text.
The current model does not expose a validated CDS phase. Emit `.` for phase in Phase 4 and
document that choice; do not infer phase from a raw `codon_start` qualifier without a dedicated
translation-aware contract.

## 5.8 Compound and cross-origin locations

A single flat interval can misrepresent a compound feature.

Recommended GFF3 behavior:

- emit one row per `LocationPart`;
- reuse a stable logical feature ID with a part index attribute, or use parent/child IDs;
- preserve strand;
- mark `cross_origin=true` when applicable.

Recommended BED behavior:

- emit one BED6 row per location part;
- add `part=1/2` style information to the name or an optional extra column in a documented BED-like extended mode.

Avoid BED12 in the first pass because circular cross-origin features do not map cleanly onto one linear block structure.

## 5.9 Localized review findings

Claude's browser-track use case is strongest when warnings such as `partial_cds_at_edge` can themselves be viewed as intervals.

Do **not** parse diagnostic message text to rediscover the affected feature.

Instead, introduce an internal richer review result:

```python
@dataclass(slots=True, frozen=True)
class ReviewFinding:
    diagnostic: Diagnostic
    location: Location | None
    entity_type: str | None
    entity_id: str | None
```

Then:

```python
def review_findings(record: Record) -> list[ReviewFinding]:
    ...


def review_record(record: Record) -> list[Diagnostic]:
    return [finding.diagnostic for finding in review_findings(record)]
```

This preserves the existing public `review_record()` return type while giving GFF3/BED a structured way to emit diagnostic tracks.

For example:

```text
review_diagnostic  partial_cds_at_edge  EBIMEA_021660  warning
```

This is preferable to adding coordinate fields directly to `Diagnostic` and changing every current JSON diagnostic entry solely for browser export.

## 5.10 GFF3/BED tests

Required tests:

- internal `[0, 100)` -> GFF3 `1..100`;
- internal `[0, 100)` -> BED `0 100`;
- negative-strand feature;
- fuzzy/partial feature;
- two-part compound feature;
- cross-origin feature;
- missing locus tag gets stable fallback ID;
- escaping of semicolon, equals, comma, spaces, and percent signs;
- localized `partial_cds_at_edge` finding appears at the CDS coordinates;
- repeated render produces byte-identical output.

---

# 6. Phase 5 - provenance and reproducibility manifest

## 6.1 Goal

Answer two distinct questions:

1. **What biological/input file did this review come from?**
2. **What antiSMASH/review software and metadata produced the annotations being compared?**

The package already stores:

- `source_path`
- `source_sha256`
- `antismash_version`
- its own parser/package version in JSON exports

The sixth pass should make this provenance explicit and extensible.

### Phase 5 entry and exit gate

Entry: a survey has captured structured-comment examples from the supported antiSMASH versions
that are actually available. Do not invent antiSMASH 5/7 fixtures or keys from memory.

Exit:

- raw keys and repeated values are preserved without false normalization;
- absent metadata remain unknown rather than equal;
- source-level manifest output works without changing `Record` JSON;
- only an explicitly approved embedding of provenance bumps the record schema;
- comparison deltas are attached at a defined level (matched record or whole input), use
  `bool | None` consistently for unknown, and are covered by JSON/Markdown tests; and
- private aggregate/region bundles remain optional and ignored.

Implement Phase 5 in two subphases: first ship a manifest exporter built from existing
per-record source/version fields plus generalized metadata parsing; only then decide whether a
typed provenance field belongs on `Record`. The manifest does not need a record-schema bump.

## 6.2 Generalize antiSMASH structured-comment parsing

Current `_antismash_version()` reads only `Version` from either Biopython's `structured_comment` or the raw `##antiSMASH-Data-START##` block.

Replace the one-field parser with:

```python
def _antismash_metadata(record: SeqRecord) -> dict[str, str]:
    ...
```

Requirements:

- prefer the structured Biopython representation when present;
- fall back to parsing the antiSMASH comment block;
- retain **all** key/value pairs;
- normalize keys only for known-field lookup;
- preserve original/raw keys and values for forward compatibility;
- do not infer missing database versions from the antiSMASH software version.

Then derive:

```python
version = metadata.get(normalized_version_key)
```

for backward compatibility.

## 6.3 Typed provenance model

Add:

```python
@dataclass(slots=True)
class AntiSmashProvenance:
    version: str | None
    run_date: str | None
    pfam_version: str | None
    detection_rule_set_version: str | None
    database_versions: dict[str, str]
    raw_fields: dict[str, tuple[str, ...]]
```

and:

```python
    Record.antismash_provenance: AntiSmashProvenance
```

Keep the existing `Record.antismash_version` for one compatibility cycle, populated from `antismash_provenance.version`, or replace it only with an explicit documented migration.

Because `Record` is part of the supported typed API, do not add a required constructor field
without a default. Prefer an immutable/empty provenance default or a separate optional field
until the schema migration is released; update direct `Record(...)` construction tests before
making the field mandatory.

Because the first implementation emits provenance through the dedicated manifest rather than
embedding it in the record JSON envelope, `dumps_records()` removes the internal typed field and
the record schema remains `0.2.0`. A future deliberate embedding must bump the record schema in
the same change that alters the serialized model.

## 6.4 Do not hard-code Pfam/rule-set keys before fixture survey

The exact metadata written by antiSMASH has changed across versions and workflows.

Before implementing aliases, collect representative structured comments from:

- an antiSMASH 5-era fixture;
- antiSMASH 7;
- antiSMASH 8.0.x bacterial output;
- if relevant, fungal output.

Build a tiny key-alias registry from observed values, for example conceptually:

```python
KNOWN_PROVENANCE_KEYS = {
    "version": (...),
    "run_date": (...),
    "pfam_version": (...),
    "detection_rule_set_version": (...),
}
```

Unknown keys remain in `raw_fields` and are never discarded.

## 6.5 Manifest exporter

Add:

```text
antismash_review/exporters/provenance.py
```

Suggested formats:

```bash
antismash-review inspect result/ --format provenance-json
antismash-review inspect result/ --format provenance-tsv
```

Suggested JSON shape:

```json
{
  "schema_name": "antismash-review-provenance",
  "schema_version": "0.1.0",
  "review_tool": {
    "name": "antismash-review",
    "version": "0.1.0"
  },
  "inputs": [
    {
      "source_path": "...",
      "source_sha256": "...",
      "records": ["contig_1"],
      "antismash": {
        "version": "8.0.4",
        "run_date": "...",
        "pfam_version": null,
        "detection_rule_set_version": null,
        "database_versions": {},
        "raw_fields": {}
      }
    }
  ]
}
```

Deduplicate run metadata by source file/hash rather than repeating identical blocks for every record in an aggregate GenBank file.

## 6.6 Comparison provenance deltas

This is where provenance becomes operationally valuable.

Extend comparison output with a result such as:

```python
@dataclass(slots=True, frozen=True)
class ProvenanceDelta:
    antismash_version_changed: bool | None
    left_antismash_version: str | None
    right_antismash_version: str | None
    pfam_version_changed: bool | None
    detection_rule_set_changed: bool | None
    differing_raw_fields: tuple[str, ...]
```

Then the Markdown comparison can say:

```text
Provenance differs between runs:
- antiSMASH: 8.0.3 -> 8.0.4
- Pfam metadata: unchanged/unknown
- detection-rule metadata: changed/unknown
```

Use tri-state logic where metadata are absent:

For matched records, attach `ProvenanceDelta` to `RecordComparison`; for an input-level
manifest comparison, emit a separate top-level delta. Do not silently collapse differing
record-level metadata from one aggregate file into one source-level boolean.

```text
True / False / Unknown
```

Do not report “unchanged” when both sides simply lack the field.

## 6.7 Provenance tests

Test:

- old raw `##antiSMASH-Data-START##` comment;
- Biopython `structured_comment` representation;
- missing metadata;
- unknown key preservation;
- wrapped/comment-formatted values;
- multiple records from one file;
- conflicting metadata across records;
- JSON manifest determinism;
- compare result where antiSMASH version changes;
- compare result where a field is absent on one side;
- no false claim of equality when both values are unknown.

---

# 7. Phase 6 - N-way cohort mode

## 7.1 Goal

Support a strain collection such as a 47-genome compendium without manually running pairwise comparisons.

### Phase 6 entry and exit gate

Entry: Phase 0 loading/fingerprint semantics are stable and the provenance manifest can identify
each member's source files and hashes.

Exit:

- root-directory and explicit-manifest input modes have one unambiguous CLI grammar;
- member names, paths, and ordering are deterministic and duplicate names fail clearly;
- product and domain matrices define whether counts are per region, per record, or per member;
- normalized matrix keys retain canonical/raw display labels in machine-readable metadata;
- invalid members fail by default with member name and path, while any skip mode is explicit and
  reported; and
- matrix TSV/JSON outputs pass schema, determinism, and input-overwrite tests without clustering.

Primary output:

```text
rows    = strains/samples
columns = antiSMASH product classes
values  = presence/absence or counts
```

Secondary output:

```text
domain-content matrix and optional domain-Jaccard clustering
```

This should be described as a **BGC product-class repertoire** or **BGC feature pangenome summary**, not a homologous BGC pangenome in the strict gene-cluster-family sense. Product labels such as `NRPS` and `T1PKS` are classes, not evidence that two regions are homologous loci.

## 7.2 New module

Create:

```text
antismash_review/cohort.py
```

Suggested models:

```python
@dataclass(slots=True)
class CohortMember:
    name: str
    input_path: Path
    records: list[Record]
    product_counts: Counter[str]
    domain_counts: Counter[str]


@dataclass(slots=True)
class CohortResult:
    root: Path
    members: list[CohortMember]
    product_columns: list[str]
    domain_columns: list[str]
    product_matrix: list[list[int]]
    domain_matrix: list[list[int]]
    domain_jaccard: list[list[float]] | None
    cluster_order: list[str] | None
    cluster_newick: str | None
```

## 7.3 Input semantics

Recommended default directory model:

```text
cohort_root/
  strain_A/
    antiSMASH outputs...
  strain_B/
    antiSMASH outputs...
  strain_C/
    antiSMASH outputs...
```

Each immediate child directory is one cohort member and is independently passed through the existing `discover()` + shared `load_review_input()` pipeline.

CLI:

```bash
antismash-review cohort cohort_root/ --format product-matrix-tsv
```

Optional explicit manifest for more complex projects:

```text
sample<TAB>path
SM-NMZ<TAB>/data/.../SM-NMZ-antismash
C14-NMZ<TAB>/data/.../C14-NMZ-antismash
```

then:

```bash
antismash-review cohort --manifest samples.tsv --format product-matrix-tsv
```

This manifest mode is worth supporting early because real 47-strain collections rarely have perfectly uniform directory layouts.

The parser should define `cohort_root` as an optional positional argument and require exactly
one of `cohort_root` or `--manifest`. The two examples above must exercise the same loading
and overwrite-protection path.

## 7.4 Matrix semantics

### Product matrix

Provide both:

```text
binary presence/absence
integer region/product counts
```

CLI:

```bash
--value binary
--value count
```

Default: `binary` for pangenome-style visualization.

Columns must be deterministic, e.g. sorted by normalized key with the canonical/raw display label retained in metadata.

### Domain matrix

Use normalized `Domain.name` from parsed `aSDomain` evidence.

Possible values:

```text
binary presence/absence   # default for similarity
count                     # useful for architecture burden/composition
```

Do not mix `PFAM_domain` counts into the same matrix unless explicitly requested. `aSDomain` and Pfam are different evidence layers.

For Phase 6, define the aggregation units before writing code: product counts are region-level
counts across all records in one cohort member; domain counts are `aSDomain` feature counts
across that member; binary values are presence in the member. Preserve record/region detail in
JSON metadata so a matrix cell cannot be mistaken for one homologous BGC.

## 7.5 Domain-content similarity - Phase 7 clustering gate

The product/domain matrices are Phase 6. Jaccard distances, average-linkage clustering, leaf
ordering, and Newick output are Phase 7 and optional. They must not add a runtime dependency or
block Phase 6. If deterministic tie-breaking cannot be demonstrated with focused fixtures,
ship matrices without clustering.

Default similarity for optional clustering:

```text
Jaccard(A, B) = |A ∩ B| / |A ∪ B|
Distance      = 1 - Jaccard
```

Use binary domain presence for the first implementation because the interpretation is transparent.

### No new runtime dependency is necessary

For ~47 samples, a pure-Python deterministic average-linkage implementation is entirely feasible.

Implement:

```text
pairwise Jaccard distance
+ average-linkage agglomeration
+ deterministic lexical tie breaking
+ leaf order
+ optional Newick export
```

Complexity at N≈47 is trivial, and this avoids introducing NumPy/SciPy solely for cohort ordering.

If a future release adds sophisticated clustering/visualization, that can become an optional `analysis` extra.

## 7.6 CLI design

Add:

```bash
antismash-review cohort ROOT [options]
```

Suggested options:

```text
--manifest PATH
--format product-matrix-tsv|domain-matrix-tsv|json
--value binary|count
--cluster-by none|domain-jaccard
--tree-output PATH
--lenient
```

Avoid reusing pairwise `--match-by`; cohort mode is not record matching.

If one member cannot be loaded, default behavior should fail clearly with the sample name. Do not silently omit a strain from a scientific matrix. If skip behavior is desired later, give it an explicit option such as `--skip-invalid-members` and report skipped samples in the output metadata.

The cohort loader must return the union of every member's discovered GenBank, JSON, and
sidecar paths so `--output` cannot overwrite any input representation.

## 7.7 Cohort JSON schema

Add:

```python
COHORT_SCHEMA_NAME = "antismash-review-cohort"
COHORT_SCHEMA_VERSION = "0.1.0"
```

Include:

- member names;
- input paths and source hashes;
- antiSMASH provenance summaries;
- product matrix;
- domain matrix;
- normalization rules;
- similarity metric;
- clustering method;
- deterministic leaf order;
- skipped/error members if an explicit skip mode is ever enabled.

## 7.8 TSV outputs

Product matrix example:

```text
sample	NRPS	T1PKS	transAT-PKS	terpene
strain_A	1	1	0	1
strain_B	1	0	1	0
strain_C	0	1	0	1
```

Optional domain matrix:

```text
sample	AMP-binding	Condensation	PKS_KS	PKS_AT	PP-binding
...
```

Optional distance matrix:

```text
sample	strain_A	strain_B	strain_C
strain_A	0	0.33	0.75
...
```

A heatmap image is intentionally **not** part of the core sixth-pass scope. Producing clean TSV/JSON and an optional Newick tree makes the output immediately usable in R, Python, iTOL, ComplexHeatmap, seaborn alternatives, or manuscript-specific plotting scripts without burdening the parser package with plotting dependencies.

## 7.9 Cohort tests

Create three to five synthetic member directories with known products/domains and test:

- deterministic sample ordering;
- explicit manifest ordering policy;
- binary product matrix;
- count product matrix;
- domain matrix;
- duplicate sample-name error;
- malformed member reports member name/path;
- aggregate-vs-region discovery still avoids duplicate representation;
- Jaccard distances with hand-calculated expected values;
- deterministic clustering under tied distances;
- valid Newick output;
- JSON schema/version;
- 47-member synthetic smoke test to prevent accidental quadratic/cubic implementation pathologies from becoming noticeable.

---

# 8. CLI surface after the sixth pass

A coherent end state could look like:

```bash
# Existing review
antismash-review inspect result/
antismash-review inspect result/ --format json
antismash-review inspect result/ --format gene-tsv
antismash-review inspect result/ --format domain-tsv
antismash-review inspect result/ --format clusterblast-tsv

# New analytical exports
antismash-review inspect result/ --format assemblyline-tsv
antismash-review inspect result/ --format assemblyline-json

# New genome-browser exports
antismash-review inspect result/ --format gff3 --output review.gff3
antismash-review inspect result/ --format bed --output review.bed

# New provenance exports
antismash-review inspect result/ --format provenance-json
antismash-review inspect result/ --format provenance-tsv

# Existing pairwise comparison, now provenance-aware
antismash-review compare run_A/ run_B/ --format markdown
antismash-review compare run_A/ run_B/ --format json

# New cohort mode
antismash-review cohort strains/ --format product-matrix-tsv
antismash-review cohort strains/ --format product-matrix-tsv --value count
antismash-review cohort strains/ --format domain-matrix-tsv --cluster-by domain-jaccard
antismash-review cohort --manifest samples.tsv --format json --cluster-by domain-jaccard
```

If the `inspect --format` choices become unwieldy, a later release can split exporters into named subcommands. For this pass, extending the established format pattern is the lowest-friction compatibility path.

---

# 9. Suggested file layout

```text
antismash_review/
  __init__.py
  _version.py
  assemblyline.py          # new: module/monomer chain reconstruction
  architecture.py          # new: domain expectation registry + assessments
  chemistry.py             # new: curated monomer masses / constants
  clustering.py            # new: optional deterministic domain-Jaccard clustering
  cli.py                   # thinner orchestration
  clusterblast.py
  cohort.py                # new: N-way cohort analysis
  compare.py               # refactored to shared fingerprints/provenance
  discovery.py
  fingerprints.py          # new: product/domain/diagnostic feature vectors
  genbank.py
  loading.py               # new: shared discover->parse->enrich loader
  locations.py             # new: shared overlap/containment helpers
  models.py
  review.py                # review_findings + diagnostic projection
  schema.py
  exporters/
    assemblyline_json.py    # new
    assemblyline_markdown.py# new
    assemblyline_table.py   # new
    bed.py                  # new
    cohort_json.py          # new
    cohort_table.py         # new
    compare_json.py
    compare_markdown.py
    entity_tables.py
    gff3.py                 # new
    json_export.py
    markdown.py
    provenance.py           # new
    tables.py
```

This keeps five concepts distinct:

```text
parsing       -> genbank.py / clusterblast.py
loading       -> loading.py / discovery.py
analysis      -> assemblyline.py / architecture.py / fingerprints.py / cohort.py
review        -> review.py
presentation  -> exporters/* / cli.py
```

That separation will matter as the tool grows.

---

# 10. Public Python API additions

The fifth pass explicitly established a typed top-level API. Extend it conservatively.

Candidates for `antismash_review.__all__`:

```python
predict_assembly_lines
assess_architecture
review_findings
render_gff3
render_bed
render_provenance_json
build_cohort
```

Do not re-export every internal registry or helper.

Add each name only at its phase's exit gate. In particular, do not expose `build_cohort` or
`render_provenance_json` before their input/error/schema contracts are tested. Public functions
must have concrete signatures in the implementation patch, not only a name in `__all__`.

Suggested stable analytical API:

```python
from pathlib import Path
from antismash_review import (
    assess_architecture,
    parse_genbank,
    predict_assembly_lines,
    review_findings,
)

record = parse_genbank(Path("region.gbk"))[0]
chains = predict_assembly_lines(record)
architecture = assess_architecture(record)
findings = review_findings(record)
```

Add import-level tests and continue shipping `py.typed`.

---

# 11. Scientific interpretation contract

Add a dedicated section to `references/semantic-contract.md` covering the new inference layer.

At minimum document the following non-claims.

## 11.1 Assembly-line prediction

- Module order is inferred from antiSMASH module annotations and genomic/protein order; cross-CDS order can be heuristic.
- `monomer_pairings` are antiSMASH predictions, not experimentally confirmed substrate incorporation.
- A-domain specificity predictors have non-zero error rates and may disagree.
- Iterative modules break the one-module-one-incorporation assumption.
- A thioesterase/release domain does not uniquely establish hydrolysis versus cyclization.

## 11.2 Mass estimation

- Reported mass is a modeled core-scaffold mass under explicit assumptions.
- Fatty-acyl tails, sugars, halogens, oxidation states, methylation, crosslinks, macrocyclization, and other tailoring can alter the final metabolite mass.
- Lipopeptides such as serrawettin-family compounds can therefore differ substantially from a peptide-only core estimate.
- Unknown monomers or hybrid PKS chemistry invalidate a complete core-mass estimate unless explicitly modeled.
- The value is a hypothesis generator for LC-HRMS prioritization, not a compound-identification criterion.

## 11.3 Architecture score

- The score is required-domain-slot coverage, not probability of biosynthetic activity.
- Missing domains can reflect contig boundaries, gene-calling issues, split genes, unusual architecture, or detection failure.
- Presence of the expected domains does not prove that the pathway is functional.

## 11.4 Cohort matrix

- Product-class presence is not BGC homology.
- Domain-content similarity is not phylogenetic distance.
- A product-class matrix is a repertoire summary, not a gene-cluster-family pangenome unless a separate homology/clustering method is applied.

These warnings should be concise in normal output but explicit in the semantic contract.

---

# 12. Backward compatibility

## Preserve

- existing `inspect` behavior;
- existing `compare` matching modes;
- existing coordinate semantics;
- existing strict/lenient parser behavior;
- current entity TSV column meanings;
- native antiSMASH JSON remaining sidecar/enrichment evidence rather than a full primary-input adapter;
- `review_record(record) -> list[Diagnostic]`.

## Changes requiring versioning

- adding `Record.antismash_provenance` changes record JSON because of `asdict()`;
- embedding provenance into `ComparisonResult` changes comparison JSON;
- adding fields to existing TSVs should be avoided unless a deliberate schema change is documented.

Prefer new formats over silently appending columns to existing stable entity TSVs.

---

# 13. Quality gates

Every phase should pass the existing gates:

Before committing this plan or any phase implementation, validate UTF-8 decoding, fenced-block
balance, and Markdown links. Documentation hygiene is part of Phase 0.

```bash
python -m ruff check .
python -m ruff format --check antismash_review tests
python -m mypy antismash_review
python -m pytest -q
python -m pytest --cov=antismash_review --cov-report=term-missing
```

On this Windows checkout, use a disposable OS temp base instead of a repository-local
pytest basetemp, and use the checkout-safe Git invocation:

```powershell
$phasePytestBase = Join-Path $env:TEMP ("antiSMASH-review-phase-" + [guid]::NewGuid())
python -m pytest -p no:cacheprovider --basetemp=$phasePytestBase -q
python -m pytest -p no:cacheprovider --basetemp=$phasePytestBase `
  --cov=antismash_review --cov-report=term-missing -q
git -c safe.directory=* diff --check
```

Add package-level verification after the new public API/exporters land:

```bash
python -m build
python -m pip install --force-reinstall dist/*.whl
python -c "import antismash_review; print(antismash_review.__version__)"
```

and a wheel-content test/check that `py.typed` remains present.

For output formats, add deterministic golden-text tests where appropriate:

- GFF3;
- BED;
- cohort TSV;
- provenance JSON;
- assembly-line TSV/JSON.

Avoid enormous fixture snapshots; use focused synthetic records whose expected biology and coordinates are obvious by inspection.

---

# 14. Suggested commit sequence

Keep commits reviewable and independently testable.

```text
0. test/fix: freeze fifth-pass baseline, warning policy, and private-fixture excludes
1. refactor: extract shared review-input loader and fingerprints
2. feat: add evidence-only NRPS/PKS assembly-line predictions
3. docs/test: freeze assembly-line evidence and ambiguity contract
4. feat: add optional curated chemistry and core-mass candidates
5. feat: add domain-aware architecture assessments and diagnostics
6. feat: add localized review findings
7. feat: add GFF3 and BED exports
8. feat: capture antiSMASH provenance metadata and manifest exports
9. feat: surface provenance deltas in compare
10. feat: add cohort product/domain matrices
11. feat: optionally add deterministic domain-Jaccard clustering
12. docs: extend semantic contract and skill workflow
13. test: add cross-version/provenance and real-output integration coverage
```

If schema changes are grouped, perform the version bump in the commit that first changes the serialized model and keep subsequent commits on that new schema version.

Do not merge commits 4-11 as one feature branch. Each commit must leave the package in a
testable state, and the optional clustering commit must be independently revertible.

---

# 15. Definition of done

The sixth pass is complete when all of the following are true.

For implementation, use these independent stopping points:

### Phase 0 gate

- [x] current `f112dd9` baseline is recorded;
- [x] private integration inputs are ignored or external;
- [x] loader and fingerprints are extracted without pairwise-output drift;
- [x] no schema version changes are required.

### Phase 1 gate

- [x] local assembly-line chains are deterministic;
- [x] monomer pairings, `X`, non-proteinogenic calls, reverse strands, and multi-CDS modules
      remain explicit;
- [x] no mass estimate is used to infer ordering, and no cross-CDS proximity heuristic is emitted;
      Phase 2 mass candidates remain separate and null for unsupported chemistry.

### Phase 2 gate

- [x] chemistry aliases have source-backed fixtures and independent mass checks;
- [x] linear and cyclic candidates are separate fields;
- [x] unresolved/hybrid chemistry returns null full-core candidates with machine-readable caveats.

### Phase 3 gate

- [x] domain slots and product rules are source-verified;
- [x] legacy diagnostics remain compatible;
- [x] assessments are scoped and exported independently of diagnostic text.

### Phase 4 gate

- [x] review findings have structured locations;
- [x] GFF3/BED `seqid`, ID, phase, escaping, compound, and rebased-region contracts are tested.

### Phase 5 gate

- [x] raw provenance fields are preserved;
- [x] manifest output works without an unnecessary `Record` schema bump;
- [x] comparison uses explicit unknown/changed/unchanged semantics.

### Phase 6 gate

- [x] product/domain matrices are deterministic and independently useful;
- [x] root-versus-manifest input semantics are unambiguous;
- [x] invalid-member behavior is explicit.

### Phase 7 gate

- [x] clustering is optional, deterministic, tested under ties, and independently revertible.

## Assembly-line interpretation

- [x] `Module.monomer_pairings` are parsed into typed substrate/monomer calls.
- [x] Reverse-strand and multi-CDS modules are handled correctly.
- [x] Ambiguous cross-CDS assembly lines remain explicitly ambiguous.
- [x] Release-domain presence is reported without assuming cyclization.
- [x] Core mass is emitted only when the modeled chemistry is complete enough.
- [x] Unresolved tails/monomers/tailoring are machine-readable caveats.

## Architecture assessment

- [x] Canonical T1PKS, transAT-PKS, and NRPS expectations are distinct.
- [x] trans-AT systems are not falsely penalized for lacking cis AT domains.
- [x] starter/final/iterative/incomplete module states are respected.
- [x] architecture score semantics are documented as parsed-domain coverage only.
- [x] unsupported product classes are `not_applicable`, not forced into a bad rule.

## Interoperability

- [x] GFF3 uses one-based inclusive coordinates.
- [x] BED retains zero-based half-open coordinates.
- [x] compound/cross-origin features do not become false contiguous intervals.
- [x] CDS/domain/module/collection tracks are deterministic.
- [x] localized review findings can be visualized in a genome browser.

## Provenance

- [x] all antiSMASH structured-comment fields are preserved.
- [x] known version/run/database fields are normalized when present.
- [x] absent metadata remain unknown rather than inferred.
- [x] source SHA-256 and review-tool version are exported.
- [x] comparison reports provenance differences/unknowns safely.

## Cohort

- [x] a directory or manifest of strains can be loaded reproducibly.
- [x] product presence/absence and count matrices are deterministic.
- [x] domain matrices are available separately.
- [x] optional domain-Jaccard clustering is deterministic.
- [x] cohort output never silently drops an invalid strain.
- [x] the documentation does not equate product-class presence with homologous BGC families.

## Repository quality

- [x] ruff passes.
- [x] ruff format check passes.
- [x] strict mypy passes.
- [x] pytest passes.
- [x] coverage does not materially regress (129 tests, 91% package coverage).
- [x] wheel installs and retains `py.typed`.
- [x] `SKILL.md` and `references/semantic-contract.md` describe all new user-visible behavior.

---

# 16. What I would deliberately leave for a seventh pass

The following are attractive, but including them now would make the sixth pass much harder to validate scientifically:

1. **Direct LC-HRMS peak-list matching** with ppm tolerances/adduct enumeration.
2. **Full PKS chemical mass reconstruction** including extender-unit selection, KR/DH/ER reduction state, cyclization, and starter-unit chemistry.
3. **Automatic fatty-acyl-tail inference** for lipopeptides.
4. **BGC homology families** across strains using sequence similarity/synteny rather than product labels.
5. **MIBiG-aware cohort clustering** or BiG-SCAPE-like similarity.
6. **Interactive plotting/heatmaps** inside the package.
7. **Native antiSMASH JSON as a complete primary input adapter**.
8. **Automatic metabolite identity claims** from predicted monomer chain + mass.

Those are valuable directions, but the sixth pass should first make the evidence model, interpretation layer, genome-browser export, provenance, and cohort summary robust enough to support them.

---

# 17. Recommended priority for actual bench/research value

If implementation time is limited, prioritize in this order:

```text
0. baseline and shared loader/fingerprints
1. evidence-only assembly-line monomer chain
2. conservative mass estimate, only after its scientific gate
3. architecture-aware completeness
4. localized findings plus GFF3/BED
5. provenance manifest
6. cohort product/domain matrices
7. optional clustering
```

For metabolite follow-up, the first three can directly help triage LC-HRMS candidates and decide whether an antiSMASH BGC is structurally plausible before investing in extraction/purification work. GFF3/BED then closes the loop with ONT read evidence, while provenance and cohort mode make the tool more reproducible and scalable across larger comparative projects.

The most important design principle throughout is simple:

> **Keep parsed antiSMASH evidence, derived inference, and biological claims as three visibly separate layers.**

That separation is what will let `antiSMASH-review` become more analytically ambitious without becoming overconfident.
