# antiSMASH-review patch implementation plan

**Repository:** `WhyAdr/antiSMASH-review`  
**Validated baseline:** `main` commit `753b223cbf565bf639a98962062dbf7782da1138` (`feat: complete sixth antiSMASH review extension`)  
**Primary purpose:** correct the Phase-2 peptide-mass semantics, harden the most visible export and real-output regression surfaces, and close a topology-semantics seam without expanding biological claims.  
**Secondary purpose:** define safe boundaries for a later LC-HRMS matcher and BGC-family work.

---

## 1. Executive verdict on the external assessment

The assessment is strong and its main scientific finding is correct. The patch should, however, be narrower and more explicit than the proposed follow-on ideas.

| Assessment item | Verdict | Refinement |
|---|---|---|
| Core-mass estimator is systematically wrong | **CONFIRMED — critical correctness bug** | The registry stores **residue** formulas/masses while the API and polymerization equation claim **free monomer** masses. For a chain of `n` supported amino acids, the current linear and cyclic candidates are each too low by exactly `n × H2O`. |
| Existing mass tests can pass while chemistry is wrong | **CONFIRMED** | `test_assemblyline.py` reproduces the same residue formulas as production code, then applies the same dehydration equation. `test_chemistry.py` tests internal consistency, not external chemical truth. |
| Markdown exporter has an important test gap | **CONFIRMED in substance** | There is no dedicated `tests/test_markdown.py`, and the ClusterBlast rendering branch needs direct tests. Re-measure the exact coverage percentage locally before quoting `43%`; this review did not independently execute coverage. Also note that `exporters/markdown.py` is the default **inspect** renderer; `compare` uses a separate `compare_markdown.py`. |
| Real biological integration coverage does not match the stated policy | **VALID concern** | `SKILL.md` explicitly requires tests using local biological files to skip cleanly when absent. Add a proper opt-in integration harness regardless of whether every historical test file contains zero `skip` calls. Do not make private biological files required for the normal suite. |
| `Location.cross_origin` vs `Record.topology` is a semantic seam | **CONFIRMED seam, not necessarily a bug** | Keep `Location.cross_origin` as **structural evidence** inferred from compound parts touching both record ends. Do **not** make it disappear merely because topology metadata are absent or wrong. Instead centralize topology interpretation and explicitly reconcile/document the two concepts. |
| Freeze all model objects after parsing | **Do not include in this patch** | `frozen=True` is not deep immutability when fields contain mutable lists/dicts. A correct immutable-model refactor requires container-type/API changes and deserves a separate design pass. |
| Direct LC-HRMS matching is high ROI | **AGREE, after the hotfix** | Implement only after mass semantics are corrected and regression-tested. Prefer a dedicated `match-hrms` analysis surface rather than pretending a peak list is merely another exporter format. |
| Extend current domain-Jaccard clustering into a lightweight “GCF” call | **REJECT the GCF label** | Current clustering is strain-level aggregate domain presence. GCFs are BGC-level homology/similarity claims. A dependency-free domain/order profile cluster can be useful, but must be named accordingly; true GCF work needs sequence-level evidence or an external method such as BiG-SCAPE. |
| Parallel cohort loading / HTML exporter | **Low priority** | Useful only after correctness and regression coverage. Neither belongs in the hotfix branch. |

### Patch ordering & execution status (v0.2.0)

```text
COMPLETED IN v0.2.0:
P0  Reproduce and lock the baseline (129 -> 144 tests passing)
P1  Correct free-monomer chemistry, centralize water mass, and guard 1-residue cyclic mass
P2  Add external-truth chemistry regression tests (1e-6 tolerance, anti-regression tests)
--  Embed provenance in record JSON envelope (schema 0.3.0)
--  Codify and test not_applicable architecture policy for unmodeled product classes
P6  Documentation, semantic contract, v0.1.0 mass regeneration notice, version bump to 0.2.0

DELIBERATELY DEFERRED TO SUBSEQUENT PASSES:
P3  Dedicated default-Markdown/ClusterBlast tests (tests/test_markdown.py)
P4  Optional real-output integration regression harness (tests/integration/)
P5  Centralized circular-topology helper (locations.py / compare.py)
P7  LC-HRMS candidate matching (separate branch)
P8  BGC-level similarity/GCF strategy (separate branch)
```

Do **not** combine P7/P8 with P1-P6.

---

# 2. P0 — baseline reproduction before modification

The current `6th-extend-plan.md` records `129` tests and `91%` package coverage, but this implementation patch must independently reproduce the baseline before editing code.

Run:

```bash
python -m pip install -e '.[test]'
python -m ruff check .
python -m ruff format --check antismash_review tests
python -m mypy antismash_review
python -m pytest -q
python -m pytest --cov=antismash_review --cov-report=term-missing -q
python -m antismash_review --help
```

Capture in the PR/agent log:

- Python version;
- Biopython version;
- exact passing/skipped test count;
- total package coverage;
- coverage for `antismash_review/exporters/markdown.py`;
- package version and assembly-line schema version;
- current output for the mass regression reproduction below.

### Mandatory bug reproduction

Before the fix, demonstrate the failure with an isolated calculation or focused test:

```python
from antismash_review.chemistry import free_monomer_mass

ser = free_monomer_mass("Ser")
leu = free_monomer_mass("Leu")
```

At the baseline, those values are residue masses approximately:

```text
Ser  87.03202840472
Leu 113.08406397853
```

but the function is named and documented as returning **free-monomer** mass. The true neutral free amino-acid masses under the package's exact-isotope constants are:

```text
Ser 105.04259308875
Leu 131.09462866256
```

For a neutral linear Ser-Leu dipeptide:

```text
expected = 218.12665706728 Da
current  = 182.10552769922 Da
error    = -36.02112936806 Da = -2 × H2O
```

For a fully supported chain of length `n`, the current error is:

```text
M_reported - M_correct = -n × M_H2O
```

where, under the package's isotope constants:

```text
M_H2O = 18.01056468403 Da
```

This is a blocking correctness defect for any downstream HRMS use.

---

# 3. P1 — fix the chemistry semantics

## 3.1 Root cause

Current `antismash_review/chemistry.py` stores formulas such as:

```text
Ser = C3H5NO2
Leu = C6H11NO
Gly = C2H3NO
```

These are amino-acid **residue formulas** (free amino acid minus H2O), not hydrated neutral free amino acids.

Yet the same module exposes:

```python
FREE_MONOMER_MASSES
free_monomer_mass(...)
```

and `assemblyline._estimate_core_mass()` calculates:

```text
M_linear = sum(M_free_monomer) - (n - 1) × H2O
M_cyclic = M_linear - H2O
```

The equation is correct **only if the registry actually contains free monomer masses**.

## 3.2 Chosen fix: preserve the written/public semantics

Use the least surprising, contract-preserving solution:

> **Keep `FREE_MONOMER_MASSES` and `free_monomer_mass()` as free-monomer APIs, and change the formula registry to true neutral free amino-acid formulas.**

Do **not** silently redefine `free_monomer_mass()` to mean residue mass.

This choice is preferable for this patch because:

1. it matches the explicit Phase-2 design in `6th-extend-plan.md`;
2. it matches existing names/docstrings;
3. it leaves the existing polymerization equation conceptually correct;
4. it avoids a broader API rename/refactor in a scientific hotfix.

A future API may additionally expose explicit residue masses if useful for peptide-MS workflows, but that should be a separate, intentionally named surface such as `RESIDUE_MASSES` / `residue_mass()`.

## 3.3 Required `chemistry.py` changes

Rename the private formula table for clarity:

```python
_FREE_AMINO_ACID_FORMULAS: dict[str, dict[str, int]] = {...}
```

Populate it with neutral free amino acids. Examples:

```text
Gly C2H5NO2
Ala C3H7NO2
Ser C3H7NO3
Leu C6H13NO2
Cys C3H7NO2S
```

All 20 proteinogenic entries must be audited, not patched only for the amino acids present in current fixtures.

Add one shared constant:

```python
WATER_MONOISOTOPIC_MASS = 2 * _ATOMIC_MASSES["H"] + _ATOMIC_MASSES["O"]
```

Then import/use that constant from `assemblyline.py` rather than duplicating atomic masses there.

Keep alias behavior unchanged:

- one-letter aliases;
- `L-` / `D-` aliases;
- exact supported proteinogenic names;
- unsupported `X`, `NH2`, `ccmal`, `D-Orn`, etc. remain unresolved.

Do not add non-proteinogenic monomers merely to make more masses non-null in this patch.

## 3.4 Required `assemblyline.py` behavior

Retain:

```python
linear = sum(free_monomer_masses) - (n - 1) * WATER_MONOISOTOPIC_MASS
cyclic = linear - WATER_MONOISOTOPIC_MASS
```

provided all existing completeness gates are satisfied.

### Recommended one-residue guard

The field is explicitly named `head_to_tail_cyclic_candidate_mass_da`. For a single canonical amino acid, subtracting one water creates a formally dehydrated value but not a defensible ordinary head-to-tail cyclic peptide candidate.

Recommended behavior:

```python
cyclic = (
    linear - WATER_MONOISOTOPIC_MASS
    if total >= 2
    else None
)
```

If this guard is implemented, document it and add a test. If maintainers deliberately retain the algebraic single-residue value, document that it is only a formal dehydration candidate. Do not leave the semantics implicit.

## 3.5 No tailoring inflation

The correction must **not** expand the chemistry scope. Preserve all current nulling rules for:

- low-confidence specificity fallback;
- unresolved monomers;
- starter/acyl uncertainty;
- iterative modules;
- PKS or hybrid chemistry;
- non-proteinogenic calls not explicitly modeled;
- downstream tailoring.

The hotfix is about making the already-supported case correct, not making the supported case larger.

---

# 4. P2 — replace self-referential tests with external-truth regression tests

## 4.1 Golden values

Add test constants that are **not calculated from the production formula registry**.

Use at least the following neutral monoisotopic values (Da), rounded only at assertion time:

```text
H2O              18.01056468403
free Gly         75.03202840472
free Ser        105.04259308875
free Leu        131.09462866256
linear Gly-Gly  132.05349212541
linear Ser-Leu  218.12665706728
cyclic Ser-Leu  200.11609238325
```

Use a documented tolerance such as:

```python
pytest.approx(expected, abs=1e-6)
```

Do not calculate `expected` from `FREE_MONOMER_MASSES`, `_FORMULAS`, or the production water constant inside the same test.

## 4.2 Update `tests/test_chemistry.py`

Retain alias/unresolved tests, but add explicit external-value assertions:

```text
test_free_glycine_mass_matches_golden_value
test_free_serine_mass_matches_golden_value
test_free_leucine_mass_matches_golden_value
test_stereo_aliases_preserve_mass
```

Optional but useful invariant:

```text
for each supported amino acid:
    free mass > residue-like current-baseline value
```

Do not make that invariant the only truth test.

## 4.3 Update `tests/test_assemblyline.py`

Replace the current test that rebuilds residue formulas with hard-coded independent golden values.

Required cases:

1. **single supported amino acid** — linear mass equals free amino-acid mass;
2. **Gly-Gly** — linear mass `132.05349212541` Da;
3. **Ser-Leu** — linear `218.12665706728` Da;
4. **Ser-Leu cyclic candidate** — `200.11609238325` Da;
5. **unknown monomer** — both full-core candidates remain null;
6. **starter module** — mass remains null;
7. **hybrid PKS segment** — mass remains null;
8. **D/L alias** — mass does not change;
9. if adopting the guard, **single-residue head-to-tail cyclic candidate is null**.

## 4.4 Add a regression against the exact old failure mode

Include a named test whose failure message makes the historical bug obvious, for example:

```text
test_linear_dipeptide_does_not_double_apply_dehydration
```

The test should fail if an implementation ever returns approximately `182.1055` Da for Ser-Leu again.

## 4.5 Source/provenance comment

In the test file, add a concise comment stating that golden masses are neutral monoisotopic molecular masses and must be independently audited against an authoritative chemical/mass reference if atomic-mass constants are ever changed.

Do not cite the production table as the source of truth.

---

# 5. P3 — dedicated tests for the default inspect Markdown renderer

Create:

```text
tests/test_markdown.py
```

The goal is not “coverage for coverage's sake”; it is to freeze user-visible behavior of the default `inspect` format.

## 5.1 ClusterBlast branch matrix

Construct a typed record carrying ClusterBlast results and cover:

### A. One ordinary hit

Assert:

- section title;
- region number;
- search-type display name;
- source format/path rendering;
- total database hits;
- ranking count;
- table header;
- numeric formatting.

### B. Six or more ranked hits

Assert:

- only first five rows render;
- explicit `Showing first 5 of N hits` line appears;
- ordering is deterministic.

### C. Missing scalar values

For `None` values in:

- `num_hits`;
- `blast_score`;
- `similarity`;

assert blank Markdown table cells rather than `None` or crashes.

### D. Cell escaping

Use accession/description text containing:

- `|`;
- newline;
- tab;
- CRLF.

Assert Markdown-safe output and no broken row structure.

### E. Sidecar source paths

Cover both:

- canonical parent directories (`clusterblast/`, `knownclusterblast/`, `subclusterblast/`);
- a source outside those canonical directory names.

### F. No ClusterBlast results

Assert no empty ClusterBlast heading is emitted.

## 5.2 CLI smoke test

Add or extend a CLI test proving that:

```bash
antismash-review inspect input.gbk
```

uses Markdown by default and returns the same renderer semantics.

Do not conflate this with `compare`, whose Markdown renderer is separate.

## 5.3 Coverage gate

After adding the tests, record the new per-file coverage for:

```text
antismash_review/exporters/markdown.py
```

Target: all meaningful branches in `_escape_cell()` and the ClusterBlast rendering block should be covered. Do not impose an arbitrary 100% project-wide gate if it would incentivize low-value tests.

---

# 6. P4 — optional real antiSMASH regression corpus without committing private biology

The repository's own `SKILL.md` already states that tests using local biological integration files must skip cleanly when those files are absent. Implement that policy explicitly.

## 6.1 Directory layout

Add only the harness to git:

```text
tests/
  integration/
    __init__.py
    test_real_antismash.py
    fingerprint.py
    README.md

tools/
  freeze_integration_fingerprints.py
```

Do **not** commit BK71-I, SM-NMZ, PF_NNT, or other private biological outputs unless redistribution is separately approved.

## 6.2 Environment-variable contract

Use one explicit variable:

```text
ANTISMASH_REVIEW_INTEGRATION_MANIFEST=/absolute/path/integration-manifest.json
```

If it is unset or points to a missing file:

```python
pytest.skip(..., allow_module_level=True)
```

or an equivalent clean skip.

No failure, download, or network access in the default test suite.

## 6.3 Local manifest format

Use JSON to avoid a new dependency. Example outside the repository:

```json
{
  "cases": [
    {
      "name": "SM-NMZ-antismash8",
      "input": "/data/.../SM-NMZ-antismash",
      "expected": "/data/.../SM-NMZ-antismash.fingerprint.json"
    }
  ]
}
```

Support any number of cases; do not hard-code sample names in test code.

## 6.4 Fingerprint content

A useful regression fingerprint should be stable across machines and should not contain full sequences or absolute source paths.

Recommended fields:

```text
record_count
record_ids
region_count
candidate_cluster_count
protocluster_count
proto_core_count
gene_count
domain_count
module_count
motif_count
raw_pfam_count
deduplicated_pfam_count
product_counter
domain_counter
diagnostic_code_counter
clusterblast_result_keys
assemblyline_prediction_count
assemblyline_mass_status_counts
antismash_versions
```

Optional per-record sections are fine if deterministically sorted.

Do not include:

- absolute local paths;
- Python object reprs;
- timestamps;
- unordered dict/list iteration;
- raw nucleotide/protein sequences.

## 6.5 Fingerprint ownership

Do not overload the existing public `fingerprints.py` merely for private test snapshots unless a genuinely reusable public abstraction emerges.

Prefer a test/integration helper first. Promote it into production only when another runtime workflow needs the same contract.

## 6.6 Freeze tool

`tools/freeze_integration_fingerprints.py` should:

1. read the local manifest;
2. load each case through the normal discovery/loading path;
3. generate deterministic fingerprint JSON;
4. write only to the `expected` path explicitly specified in the local manifest;
5. never mutate repository source files automatically.

Call this an explicit **bless/freeze** action in documentation. Tests must never silently rewrite expected values.

## 6.7 What the regression should catch

The real-data test should fail on unexpected changes to:

- antiSMASH qualifier adaptation;
- region/module/domain counts;
- product/domain normalization;
- sidecar attachment;
- diagnostic production;
- assembly-line interpretation;
- corrected mass availability/value fingerprints where supported.

This is specifically meant to catch “synthetic fixture encoded our own mistaken schema assumption” failures.

---

# 7. P5 — clarify `cross_origin` versus circular topology

## 7.1 Current situation

`genbank._location()` currently sets `Location.cross_origin=True` structurally when a compound location has at least one part beginning at `0` and at least one part ending at `record_length`.

Separately, `compare.intergenic_summary()` checks:

```python
record.topology.casefold() == "circular"
```

before adding a wrap-around intergenic gap.

These answer different questions:

```text
Location.cross_origin:
    Does this feature representation span both ends of this coordinate record?

Record.topology == circular:
    Does record metadata declare the molecule circular?
```

They should not be collapsed into one boolean.

## 7.2 Design decision

**Keep `Location.cross_origin` topology-agnostic.**

Why:

- malformed or incomplete topology metadata should not erase structural feature evidence;
- extracted/rebased antiSMASH records can have unusual boundary representations;
- parser models should preserve source structure before higher-level reconciliation.

Do **not** change `_location()` to require `record.topology == "circular"` before setting `cross_origin`.

## 7.3 Centralize topology interpretation

Add a small helper in `locations.py` (or another single appropriate module):

```python
def is_circular_topology(topology: str | None) -> bool:
    return topology is not None and topology.strip().casefold() == "circular"
```

Use it from `compare.intergenic_summary()` and any future topology-aware logic.

Add tests for:

- `"circular"`, mixed case, whitespace;
- `"linear"`;
- `None`.

## 7.4 Document the two-evidence model

Update `references/semantic-contract.md`:

```text
cross_origin is inferred from feature-part structure and does not itself prove that the source molecule is circular. Record topology is separate metadata. Topology-aware analyses must reconcile both explicitly.
```

## 7.5 Optional diagnostic — only after real-data survey

Do **not** immediately create noisy warnings without checking real antiSMASH outputs.

Survey the integration corpus for:

```text
cross_origin=True + topology=linear
cross_origin=True + topology=None
```

If explicit `linear` conflicts are real and rare enough to be actionable, consider a new diagnostic such as:

```text
cross_origin_topology_conflict
```

with a non-claim that the annotation may represent an extracted/rebased record or inconsistent metadata.

Do not warn merely because topology is absent unless evidence shows that such a notice is useful.

---

# 8. P6 — compatibility, versioning, and documentation

## 8.1 Package version

The corrected mass values are user-visible scientific output. Bump the package patch version:

```text
0.1.0 -> 0.1.1
```

Update both the packaging metadata and `_version.py` through the repository's existing versioning convention.

## 8.2 Assembly-line JSON schema

The hotfix should **not change the field structure** of `MassEstimate` unless the optional one-residue guard requires only a value change (`number -> null`) within the existing nullable field.

Recommended policy for this patch:

- keep `ASSEMBLYLINE_SCHEMA_VERSION = 0.2.0` because the serialized structure is unchanged;
- rely on `parser_version = 0.1.1` to distinguish corrected outputs;
- add a release/semantic-contract note that assembly-line core masses produced by parser version `0.1.0` are invalid for supported proteinogenic chains because of the residue/free-monomer mismatch.

If the project explicitly treats semantic corrections as schema-version events, use a patch-level schema bump (`0.2.1`) consistently. Do not make an ad-hoc minor schema bump solely because numbers were wrong.

## 8.3 Semantic contract update

Under Assembly-line evidence, state explicitly:

```text
The proteinogenic registry stores neutral free amino-acid monoisotopic masses.
Linear peptide candidates subtract (n-1) waters from the free-monomer sum.
Head-to-tail cyclic candidates subtract one additional water when such a candidate is modeled.
```

Retain the existing non-claims around:

- observed metabolite identity;
- lipopeptide tails;
- tailoring;
- unresolved topology;
- hybrid PKS chemistry.

## 8.4 SKILL.md maintenance section

Add the optional integration command, for example:

```bash
ANTISMASH_REVIEW_INTEGRATION_MANIFEST=/path/to/integration-manifest.json \
python -m pytest tests/integration -q
```

Explicitly say the default suite must pass with this variable unset and the integration module skipped.

## 8.5 Release note

Add a concise note wherever release changes are tracked:

```text
Fixed a Phase-2 peptide core-mass bug in which residue formulas were mislabeled as free-monomer formulas and then dehydrated again. Proteinogenic linear/cyclic candidates from parser v0.1.0 may be lower than correct values by n×18.01056468403 Da for n incorporated supported amino acids and should be regenerated.
```

That warning is important for any already-exported LC-HRMS candidate tables.

---

# 9. Required test matrix after P1-P6

At minimum, the final suite should contain all of the following.

## Chemistry

- [ ] 20 free amino-acid formulas audited.
- [ ] Gly free monomer golden mass.
- [ ] Ser free monomer golden mass.
- [ ] Leu free monomer golden mass.
- [ ] D/L aliases mass-identical.
- [ ] unsupported monomers remain unresolved.
- [ ] water mass has a single owner.

## Assembly-line mass

- [ ] single supported monomer linear mass correct.
- [ ] Gly-Gly linear mass correct.
- [ ] Ser-Leu linear mass correct.
- [ ] Ser-Leu cyclic candidate correct.
- [ ] old double-dehydration value explicitly rejected.
- [ ] unresolved monomer -> null full-core mass.
- [ ] low-confidence specificity -> null.
- [ ] starter/acyl uncertainty -> null.
- [ ] iterative module -> null.
- [ ] PKS/hybrid -> null.
- [ ] optional one-residue cyclic guard covered.

## Markdown

- [ ] ClusterBlast ordinary rendering.
- [ ] >5 ranking truncation.
- [ ] missing score/similarity fields.
- [ ] pipes/newlines/tabs escaped.
- [ ] canonical/non-canonical source path display.
- [ ] no empty ClusterBlast heading when absent.
- [ ] default inspect CLI smoke test.

## Integration

- [ ] clean module-level skip with env var absent.
- [ ] one or more real antiSMASH bundles pass when local manifest is provided.
- [ ] fingerprints deterministic across two runs.
- [ ] expected snapshots are never auto-updated by tests.

## Topology

- [ ] structural cross-origin detection remains independent of topology.
- [ ] circular topology helper canonicalizes case/whitespace.
- [ ] intergenic wrap gap still depends on declared circular topology.
- [ ] semantic contract states that the two signals are distinct.

---

# 10. Suggested commit sequence for the Gemini agents

Keep these commits independent and reviewable:

```text
1. test: reproduce peptide core-mass regression with external golden values
2. fix: correct free-monomer formulas and centralize water mass
3. test: harden assembly-line chemistry edge cases
4. test: cover default Markdown ClusterBlast rendering
5. test: add opt-in real antiSMASH integration fingerprint harness
6. refactor/docs: centralize circular-topology handling and document cross-origin semantics
7. chore: bump patch version and document invalid v0.1.0 mass outputs
```

Do not squash the first failing regression test into the fix during agent work if preserving the red/green history is practical. It makes the scientific correction auditable.

---

# 11. Agent execution rules

Each implementing agent should follow these rules:

1. **Read the current code, not only this plan.** Confirm HEAD before editing.
2. **Do not broaden biological inference** while fixing the mass bug.
3. **Do not add chemistry aliases** without an explicit formula/reference and test.
4. **Do not derive expected test masses from production constants/tables.**
5. **Do not commit private biological integration files.**
6. **Do not rename domain-Jaccard clustering to GCF clustering.**
7. **Do not make `cross_origin` topology-dependent in the parser.**
8. **Do not “freeze” dataclasses superficially while nested containers remain mutable.**
9. Run Ruff, format check, strict mypy, pytest, coverage, CLI help, and wheel/import checks after the patch.
10. Report any output-schema change before implementing it; this patch is intended to avoid structural schema changes.

---

# 12. Post-patch extension P7 — LC-HRMS matching (separate branch)

This is the highest-value follow-on feature, but only after P1-P6 are green.

## 12.1 Architectural boundary

Do not call observed m/z support a metabolite identification.

Use terminology such as:

```text
candidate mass match
peak support
ppm-compatible candidate
```

Never:

```text
identified compound
confirmed metabolite
```

without external experimental evidence.

## 12.2 Prefer a dedicated analysis command

A peak list is an additional input, not merely an output format. Prefer:

```bash
antismash-review match-hrms result/ peaks.tsv \
  --ppm 5 \
  --polarity positive \
  --format tsv
```

over:

```text
inspect --format hrms-match-tsv
```

unless the latter also introduces an explicit mandatory `--peak-list` argument.

## 12.3 V1 peak-list contract

Keep dependencies minimal. Accept TSV/CSV first rather than native mzML.

Minimum columns:

```text
mz
```

Optional:

```text
intensity
retention_time
peak_id
```

Native mzML can be a later optional adapter if a justified dependency is accepted.

## 12.4 Candidate sources

Only match non-null corrected `MassEstimate` candidates:

```text
linear_core_mass_da
head_to_tail_cyclic_candidate_mass_da
```

Carry forward:

- record ID;
- region number;
- assembly-line/module identity;
- chain display;
- chemistry scope/caveats;
- topology candidate type.

Do not construct a mass for unresolved chains merely to increase hit count.

## 12.5 Adduct model

Start with singly charged common adducts only:

Positive:

```text
[M+H]+
[M+Na]+
[M+NH4]+
```

Negative:

```text
[M-H]-
```

**Important:** adduct ion masses must use appropriate ionic/proton mass constants. Do not reuse the neutral hydrogen-atom mass from the molecular-formula calculator for `[M+H]+`.

Source and test proton/cation masses independently before implementation.

Defer:

- multiply charged ions;
- dimers `[2M+H]+`;
- in-source fragments;
- water/ammonia losses;
- isotope-envelope scoring;
- retention-time prediction.

## 12.6 PPM calculation

Define the signed error once:

```text
delta_da  = observed_mz - theoretical_mz
delta_ppm = delta_da / theoretical_mz × 1e6
match     = abs(delta_ppm) <= ppm_tolerance
```

Sort deterministically by:

```text
abs(delta_ppm), observed_mz, candidate identity, adduct
```

## 12.7 Output columns

Recommended TSV:

```text
record_id
region_number
assembly_line_index
candidate_type
neutral_mass_da
adduct
charge
theoretical_mz
observed_mz
delta_da
delta_ppm
intensity
retention_time
chemistry_scope
caveats
```

This is candidate support only.

---

# 13. Post-patch extension P8 — BGC family work: reframe before coding

The current `clustering.py` performs **member-level** clustering from aggregate adapted-domain presence across an entire cohort member/strain.

That is useful for repertoire similarity, but it is not a BGC-family engine.

## 13.1 Do not call this a GCF

If you extend the dependency-free implementation using:

- per-BGC domain composition;
- domain order;
- coarse synteny;

name the result something like:

```text
BGC architecture similarity cluster
BGC profile cluster
```

and explicitly state that it is not a homologous GCF call.

## 13.2 A true GCF workflow needs a different unit of analysis

First construct one typed BGC entity/fingerprint per region/protocluster, not one aggregate vector per strain.

A defensible GCF-like method would additionally need sequence-level evidence, e.g.:

- homologous protein/domain sequence similarity;
- shared gene-family content;
- adjacency/synteny;
- calibrated distance/threshold behavior;
- benchmarking against known related/unrelated BGCs.

At that point either:

1. integrate/export to BiG-SCAPE and ingest its family calls; or
2. deliberately implement and benchmark a new BGC-distance method.

Do not blur the current semantic contract merely because `average_linkage_domain_clustering()` already exists.

---

# 14. Explicitly deferred items

Do not include in this patch:

- deep immutable `Record`/`Gene`/`Module` redesign;
- cohort parallelism;
- HTML exporter;
- native mzML parsing;
- full PKS mass reconstruction;
- fatty-acyl-tail inference;
- automatic metabolite identity;
- a home-grown feature called “GCF” based only on domain Jaccard.

These are reasonable future projects, but they increase scope without improving confidence in the critical chemistry fix.

---

# 15. Final definition of done

The patch is complete only when:

- [ ] baseline checks were reproduced and recorded;
- [ ] proteinogenic registry values truly represent neutral free amino acids;
- [ ] water mass has one implementation owner;
- [ ] Ser-Leu and Gly-Gly golden regressions pass;
- [ ] the historical double-dehydration failure is impossible without a test failure;
- [ ] unsupported chemistry still produces null full-core candidates;
- [ ] default inspect Markdown has dedicated ClusterBlast branch tests;
- [ ] optional private real-data integration tests skip cleanly when absent and pass when configured;
- [ ] `cross_origin` and circular topology are documented as separate evidence layers;
- [ ] package version/release note warns that v0.1.0 mass outputs must be regenerated;
- [ ] Ruff, formatting, strict mypy, pytest, coverage, CLI help, build, wheel install, and `py.typed` checks are green;
- [ ] no private biological data were committed;
- [ ] no new biological identity claim was introduced.

