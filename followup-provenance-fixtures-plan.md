# Follow-up Provenance & Fixture Hardening Plan

**Repository:** `WhyAdr/antiSMASH-review`  
**Reviewed remote:** `main`  
**Reviewed HEAD:** `fb546caa3d6cc32b9a6f4e829a8f8556c4b6a44f`  
**Current package version:** `0.3.0`  
**Current record schema:** `antismash-review / 0.3.0`  
**Review date:** 2026-08-17

---

## 1. Executive assessment

The latest patch series closes almost all of the concrete correctness issues from the previous review:

- module-level monomer interpretation is now separated from raw `/monomer_pairings`;
- malformed/unresolved pairing identities no longer collapse too aggressively;
- incomplete modules and unresolved domain references block core-mass estimates;
- architecture scoring is restricted to `nrps_pks_domains`;
- NRPS module scoring now requires explicit NRPS module typing;
- ClusterBlast merge/attachment behavior is transactional and lenient-mode failures are independently retained;
- v6/v7/v8 minimized ClusterBlast fixtures now use the correct outer module schema for the releases they claim to represent;
- multi-record attachment diagnostics have better record provenance;
- file-level parse diagnostics no longer falsely claim the first record ID;
- the semantic contract was corrected;
- CI now runs Python 3.10–3.12, Ruff, strict mypy, pytest with a 90% coverage floor, and a packaging smoke step.

The current CI run for HEAD `fb546caa...` completed successfully for Python 3.10, 3.11, and 3.12.

The remaining work is therefore narrower. I found **one real compatibility bug**, several **provenance/fixture mismatches**, and a few **regression/packaging hardening opportunities**. None of these invalidate the new assembly-line chemistry work.

The highest-priority new finding is:

> **antiSMASH 7.0.x ClusterBlast JSON uses `GeneralResults.schema_version = 2`, but `antiSMASH-review` currently rejects result schema 2.**

So the present documentation claim of supporting “antiSMASH 7.x” is too broad: 7.1 is modeled, but 7.0.x is not.

---

# 2. Source-audited ClusterBlast version matrix

The fixture/support matrix should be grounded directly in upstream serializer versions rather than inferred from major antiSMASH versions.

| antiSMASH generation | `ClusterBlastResults` / module container schema | `GeneralResults` result schema | Important change |
|---|---:|---:|---|
| 5.0–5.2 | 2 | 1 | baseline modern JSON shape |
| 6.x | 2 | 1 | same broad ranking/pairing representation |
| 7.0.x | 2 | **2** | adds optional `data_version`; no schema-3 similarity yet |
| 7.1.x | 2 | **3** | JSON similarity added |
| 8.x | 2 | **5** | later ClusterBlast representation changes; `Subject.full_name` present |
| synthetic/unknown legacy compatibility | 1 | 1 | no released-version provenance established yet |

Relevant upstream source anchors:

- `antismash/antismash@5-0-0` — `antismash/modules/clusterblast/results.py`
- `antismash/antismash@5-2-0` — `antismash/modules/clusterblast/results.py`
- `antismash/antismash@6-1-1` — `antismash/modules/clusterblast/results.py`
- `antismash/antismash@7-0-0` and `7-0-1` — `antismash/modules/clusterblast/results.py`
- `antismash/antismash@7-1-0` — `antismash/modules/clusterblast/results.py`
- `antismash/antismash@8-0-4` — `antismash/modules/clusterblast/results.py`

This matrix should become the canonical support table in tests and documentation.

---

# 3. P0 — Fix the antiSMASH 7.0.x schema-2 compatibility hole

## Problem

Current parser validation accepts:

```python
result_schema_version in {1, 3, 5}
```

and current negative tests deliberately reject `2`.

That is historically incorrect for antiSMASH 7.0.x.

Upstream 7.0.0/7.0.1 define approximately:

```python
class ClusterBlastResults(ModuleResults):
    schema_version = 2

class GeneralResults(ModuleResults):
    schema_version = 2
```

Schema 2 mainly introduced optional `data_version`.

Therefore a valid antiSMASH 7.0.x JSON result can currently be rejected as unsupported.

## Required code changes

### `antismash_review/clusterblast.py`

Change accepted result schemas from:

```python
{1, 3, 5}
```

to:

```python
{1, 2, 3, 5}
```

Do **not** flatten all versions conceptually into one undocumented allow-list. Add a comment or helper mapping:

```python
_CLUSTERBLAST_RESULT_SCHEMA_GENERATIONS = {
    1: "antiSMASH 5.x/6.x",
    2: "antiSMASH 7.0.x",
    3: "antiSMASH 7.1.x",
    5: "antiSMASH 8.x",
}
```

The mapping is documentation/validation metadata, not a requirement to branch parsing unnecessarily if the normalized fields share the same representation.

### `tests/test_clusterblast_validation.py`

Current negative parametrization that includes schema `2` must no longer reject it. Keep negative coverage for unsupported values such as `0`, `4`, and `6`, and add an explicit positive schema-2 test.

## Add a real version-specific fixture

Add:

```text
tests/fixtures/clusterblast/
    clusterblast_v7_0_schema2.json
```

Rename the current 7.x fixture to:

```text
clusterblast_v7_1_schema3.json
```

so the filename no longer implies that schema 3 represents every 7.x release.

Expected assertions:

```python
v70.module_schema_version == 2
v70.result_schema_version == 2
v70.rankings[0].similarity is None
```

For 7.1:

```python
v71.module_schema_version == 2
v71.result_schema_version == 3
v71.rankings[0].similarity == EXPECTED
```

## Documentation changes

Update:

- `references/semantic-contract.md`
- `SKILL.md` if the compatibility matrix is stated there
- `tests/fixtures/clusterblast/README.md`

Use:

```text
result schema 1 → antiSMASH 5.x/6.x
result schema 2 → antiSMASH 7.0.x
result schema 3 → antiSMASH 7.1.x
result schema 5 → antiSMASH 8.x
```

Avoid saying merely “schema 3 = antiSMASH 7.x”.

## Acceptance criteria

- A serializer-faithful antiSMASH 7.0.x JSON fixture parses.
- Existing v6, v7.1, v8 tests remain unchanged semantically.
- Unsupported schemas still fail strictly.
- CI remains green across Python 3.10–3.12.

---

# 4. P1 — Preserve ClusterBlast database/data provenance

## Problem

Upstream antiSMASH 7.0+ `GeneralResults` carries:

```python
data_version
```

and includes it in JSON when available.

Current `antiSMASH-review` drops it.

This is not just decorative metadata. It becomes scientifically important as soon as the tool derives:

- KnownClusterBlast query coverage;
- unmatched-query fraction;
- top-N KnownClusterBlast coverage;
- reference-hit repertoire Jaccard;
- comparisons between antiSMASH runs;
- any “novelty-like” triage metric.

Two runs against different ClusterBlast/MIBiG datasets should not silently appear equivalent.

## Model change

Extend:

```python
@dataclass(slots=True)
class ClusterBlastResult:
    ...
    module_schema_version: int | None = None
    result_schema_version: int | None = None
    data_version: str | None = None
```

Prefer `data_version` rather than a generic `database_version` because that mirrors the upstream field exactly.

If later KnownClusterBlast exposes richer database provenance, add separate typed fields rather than overloading this value.

## Parser changes

For each `general`, `knowncluster`, and `subcluster` section:

```python
data_version = _optional_string(section.get("data_version"), ...)
```

Preserve `None` if absent.

Also validate versioned fields that upstream serializers emit when present:

```text
section.record_id
section.search_type
```

Recommended policy:

- if absent in a minimized/legacy document, tolerate absence where historically appropriate;
- if present, require correct type;
- if `record_id` is present, require it to match the containing record/module ID;
- if `search_type` is present, require it to agree with the section key.

This provides stronger source-shape validation without making unused fields part of biological inference.

## Schema migration

`ClusterBlastResult` is nested in `Record`, and record JSON is produced with `asdict(record)`. Therefore adding `data_version` changes the public record JSON shape.

Do **not** add it silently while keeping:

```text
RECORD_SCHEMA_VERSION = 0.3.0
```

Recommended release boundary:

```text
package version:       0.4.0
record schema:         0.4.0
assemblyline schema:   0.3.0   # unchanged if assembly-line shape does not change
comparison schema:     bump only if its public shape changes
```

This is a good place for a minor package release because the change improves provenance and compatibility rather than correcting the v0.3 mass semantics.

## Export changes

Because record JSON will naturally expose `data_version`, additionally consider exposing it in:

```text
clusterblast-tsv
```

Useful provenance columns:

```text
source_format
module_schema_version
result_schema_version
data_version
```

If those provenance columns already exist partly, append rather than reorder unless a schema/versioned table contract exists.

## Tests

Add:

1. schema-2/3/5 fixture with `data_version`;
2. absence remains `None`;
3. non-string `data_version` rejected;
4. mismatching section `record_id` rejected when present;
5. mismatching section `search_type` rejected when present;
6. record JSON emits the preserved value;
7. schema version assertion updates to record `0.4.0`.

---

# 5. P2 — Split “minimal parser fixtures” from “serializer-faithful golden fixtures”

## Problem

The new fixtures are much better than changing `schema_version` inside one in-memory object, but they still serve two different roles at once.

For example the current v8 fixture is intentionally compact. It omits much of what upstream `GeneralResults.to_json()` / `RegionResult.jsonify()` emits, including fields such as:

- section `record_id`;
- `search_type`;
- `proteins`;
- optional `data_version`;
- several `ReferenceCluster` fields;
- v8 `Subject.full_name`.

Yet the README calls the files “version-faithful minimized fixtures” and specifically notes v8 `Subject.full_name`.

That is too strong for the current artifact.

The parser does not need every unused upstream field in every unit fixture, but the test suite also needs a second class of fixture that protects real serializer compatibility.

## Proposed layout

```text
tests/fixtures/clusterblast/
├── README.md
├── text/
│   └── ...
├── minimal/
│   ├── schema1_minimal.json
│   ├── schema2_minimal.json
│   ├── schema3_minimal.json
│   └── schema5_minimal.json
└── golden/
    ├── antismash_6_1_1_clusterblast.json
    ├── antismash_7_0_1_clusterblast.json
    ├── antismash_7_1_0_clusterblast.json
    └── antismash_8_0_4_clusterblast.json
```

### Minimal fixtures

Purpose:

- test field validation;
- test optional/malformed inputs;
- test normalization;
- stay small and easy to edit.

They may omit upstream fields the parser intentionally ignores.

Call them explicitly:

```text
minimal parser-contract fixtures
```

### Golden fixtures

Purpose:

- protect version compatibility against real upstream serialization shape;
- catch changes in nested subject/reference/query structures;
- document exact historical provenance.

These should be either:

**A. extracted from actual antiSMASH result JSON**, or

**B. source-reconstructed from the exact upstream serializer**, preserving all fields emitted by the relevant serializer even if values are synthetic.

If option B is used, call them:

```text
serializer-faithful reconstructed fixtures
```

not “authentic output”.

Only use “authentic” for files actually emitted by antiSMASH.

## Fixture metadata README

For every golden fixture record:

| Field | Example |
|---|---|
| antiSMASH tag | `8-0-4` |
| upstream commit SHA | exact SHA |
| serializer module | `antismash/modules/clusterblast/results.py` |
| serializer functions | `GeneralResults.to_json`, `RegionResult.jsonify` |
| related data structure | `data_structures.py:Subject` |
| source type | real output / reconstructed |
| transformations | identifiers anonymized? ranking truncated? |
| fixture SHA-256 | stable hash |
| expected normalized result | test name/reference |

For actual real outputs, prefer deterministic redaction/minimization scripts rather than hand-editing.

## Parity tests

Do not assert only that all versions happen to contain the same manually chosen subset.

Instead assert stable normalized fields separately from version-dependent fields.

### Stable normalized fields

```text
record_id
region_number
search_type
ranking accession
cluster type
hit counts
BLAST score
pairing query/subject
identity
coverage
e-value
```

### Version-dependent fields

```text
schema 1: similarity absent/None
schema 2: data_version may be present; similarity absent/None
schema 3: similarity serialized
schema 5: similarity serialized; newer Subject representation tolerated
```

This makes the tests document historical differences rather than erase them.

---

# 6. P3 — Resolve diagnostic ownership instead of using record 0 as a sink

There are two related provenance levels:

```text
input/source-level diagnostics
record-level diagnostics
```

They should not be conflated.

## P3A — Duplicate-result diagnostics are still stored on record 0

The final patch correctly improved `clusterblast_attach_failed` routing.

However, current loading does:

```python
diag_sink = records[0].diagnostics if (lenient and records) else None
merged = merge_clusterblast_results(
    ...,
    diagnostics=diag_sink,
)
```

`merge_clusterblast_results()` creates a diagnostic with:

```python
record_id=result.record_id
```

but the object is still physically appended to `records[0].diagnostics`.

Example:

```text
Record A
Record B

duplicate ClusterBlast result for B
```

can produce:

```text
Record A.diagnostics:
    code = clusterblast_duplicate_result
    record_id = B
```

This is internally contradictory.

It matters downstream because:

- record JSON embeds `record.diagnostics`;
- `review_record(record)` begins with the record’s stored diagnostics;
- pairwise comparison fingerprints diagnostics **per matched record**.

Thus a duplicate involving B can be counted as a change in A.

## Near-term fix

Collect merge diagnostics separately:

```python
merge_diagnostics: list[Diagnostic] = []

merged = merge_clusterblast_results(
    text_results,
    json_results,
    lenient=lenient,
    diagnostics=merge_diagnostics,
)
```

Then route each diagnostic:

```python
matches = [
    record
    for record in records
    if diagnostic.record_id in {record.record_id, record.name}
]
```

Policy:

- exactly one match → append to that record;
- zero matches → input-level diagnostic;
- multiple matches → input-level diagnostic because attribution is ambiguous.

Never silently fall back to record 0 for a diagnostic whose owner is another record.

## Required regression

Construct:

```text
Record A
Record B
duplicate ClusterBlast result for B
```

and assert:

```python
not any(d.code == "clusterblast_duplicate_result" for d in A.diagnostics)
any(d.code == "clusterblast_duplicate_result" for d in B.diagnostics)
```

Also compare A/B before/after and prove B, not A, gets the diagnostic delta.

---

# 7. P3B — File-level parse failures need a first-class input diagnostic container

The final patch appropriately sets:

```python
record_id=None
```

for sidecar parse failures whose record identity is unknown.

But they are still physically appended to:

```python
records[0].diagnostics
```

That is better than falsely labeling record 0, but the ownership problem remains.

Example:

```text
multi-record GenBank
+
corrupt JSON sidecar that cannot be parsed far enough to identify a record
```

becomes a diagnostic nested under record A even though it belongs to the input/result directory as a whole.

This can pollute record-level comparisons.

## Preferred architecture

Extend:

```python
@dataclass(slots=True)
class LoadedReviewInput:
    root: Path
    records: list[Record]
    input_paths: set[Path]
    input_diagnostics: list[Diagnostic]
```

Use `input_diagnostics` for:

- malformed ClusterBlast sidecar before record identity is known;
- duplicate result when record ownership cannot be resolved uniquely;
- attachment failure when no record ID/name matches;
- other future directory-level enrichment failures.

Keep `Record.diagnostics` exclusively for record-scoped diagnostics.

## Export semantics

Top-level JSON already has a document-level `diagnostics` collection.

Refactor it deliberately:

```text
records[].diagnostics
    = record-scoped source diagnostics only

top-level diagnostics
    = input-level diagnostics
      + all record-level diagnostics
      + derived review diagnostics
```

Avoid duplicating the same diagnostic twice unless the schema explicitly says the top-level collection is an aggregate index.

If top-level diagnostics remain an aggregate index, document that clearly.

## Comparison semantics

Input-level diagnostics should **not** be inserted into `RecordComparison.new_diagnostics` / `resolved_diagnostics`.

Options:

1. add separate fields to `ComparisonResult` such as `left_input_diagnostics` and `right_input_diagnostics`, and compare them separately; or
2. retain them as metadata but do not calculate deltas yet.

If the comparison JSON public shape changes, bump `COMPARISON_SCHEMA_VERSION`. Do not hide a public-shape change under `0.2.0`.

## Compatibility wrapper

`cli.load_review_records()` currently returns `(records, input_paths)`.

It may stay as a compatibility wrapper and intentionally discard input diagnostics for old Python callers, but new internal CLI paths should use `LoadedReviewInput`.

Document that limitation until the old wrapper is deprecated.

---

# 8. P4 — Reassess module-container schema 1 support

Current parser intentionally accepts outer ClusterBlast module schemas:

```python
{1, 2}
```

The new fixture calls schema 1:

```text
Legacy/custom
```

That is appropriately more cautious than tying it to antiSMASH 6 or 7.

However, the current review found:

- antiSMASH 5.0.0: outer schema 2;
- antiSMASH 5.2.0: outer schema 2;
- antiSMASH 6.x: outer schema 2;
- antiSMASH 7.x: outer schema 2;
- antiSMASH 8.x: outer schema 2;
- an inspected pre-v5 2018 ClusterBlast implementation also already used outer schema 2.

So schema-1 outer-container support currently lacks a known released-version provenance anchor.

## Recommended action

Do **not** remove support immediately; permissive compatibility is harmless if parsing semantics are safe.

Rename/document it as:

```text
synthetic compatibility fixture — no upstream release provenance established
```

For example:

```text
clusterblast_compat_module_schema1.json
```

README wording:

> This fixture exercises an intentionally tolerated module-container schema 1 shape. No released antiSMASH tag has yet been identified as its provenance; it must not be used as evidence that a specific antiSMASH release emitted this schema.

## Optional history audit

Search older tags/branches/history for the exact commit that introduced `ClusterBlastResults.schema_version = 2`.

If a genuine schema-1 producer is identified:

- pin tag/commit;
- reconstruct/extract a golden fixture;
- update README.

If none is found after a bounded audit, consider whether schema-1 support is worth retaining indefinitely.

---

# 9. P5 — Make the wheel smoke test genuinely isolated

## Current state

CI is now meaningful and is passing on Python 3.10–3.12.

The packaging smoke step currently runs in the same environment where the project has already been installed editable:

```bash
pip install -e ".[test]"
...
python -m build
pip install dist/*.whl
antismash-review --help
```

This is better than having no packaging test, but it does not conclusively prove that the built wheel is independently installable.

Potential masking mechanisms:

- the same package/version is already installed editable;
- pip may regard the wheel as already satisfied;
- repository current working directory can shadow installed package imports.

## Preferred smoke test

Use a clean venv and run from outside the repository:

```bash
python -m build

python -m venv /tmp/antismash-review-wheeltest
/tmp/antismash-review-wheeltest/bin/python -m pip install --upgrade pip
/tmp/antismash-review-wheeltest/bin/python -m pip install dist/*.whl

cd /tmp
/tmp/antismash-review-wheeltest/bin/antismash-review --help
/tmp/antismash-review-wheeltest/bin/python - <<'PY'
import antismash_review
from pathlib import Path

print(antismash_review.__version__)
print(Path(antismash_review.__file__).resolve())
PY
```

Assert the imported path is inside the temporary venv/site-packages and not the GitHub workspace.

An even cleaner option is a separate `package` job depending on the test matrix.

## Optional addition

Run:

```bash
python -m pip check
```

inside the wheel-test venv.

---

# 10. P6 — Add the deferred real-data integration harness

The synthetic suite is now strong enough that the next confidence gain comes from controlled real outputs.

No opt-in integration harness was found at current HEAD.

Recommended environment contract:

```text
ANTISMASH_REVIEW_INTEGRATION_MANIFEST=/path/to/manifest.tsv
```

Example manifest:

```text
name    path
sm_zpg19    /private/results/SM-ZPG19
bk71_i      /private/results/BK71-I
```

Tests should skip cleanly when the environment variable is absent.

## Freeze only deterministic fingerprints

Do not commit private sequences or full result directories.

For each integration case, freeze:

```text
record IDs
antiSMASH version
source-sidecar basename/type
ClusterBlast schema/data_version
region/product counts
domain/module counts
assembly-line module count
raw pairing counts
interpreted incorporation counts
pairing_status
mass-null/non-null state
diagnostic codes
selected KnownClusterBlast accessions where licensing/privacy permits
```

Avoid:

- absolute paths;
- timestamps;
- raw nucleotide/protein sequences;
- user-specific filesystem details.

## High-value real regression

The existing cross-CDS duplicated `Orn -> D-Orn` case should be represented as a real-data integration fingerprint:

```text
raw monomer pairings = 2
interpreted incorporation slots = 1
pairing_status = identical_duplicate
cross_cds_duplicate_monomer_pairing present
```

That locks the scientific motivation for the v0.3 semantic change.

---

# 11. P7 — Add dedicated Markdown renderer regression tests

A repository search still found no dedicated:

```text
tests/test_markdown.py
```

Existing CLI/export tests provide some indirect Markdown coverage, but renderer regressions are easy to miss when table shapes and caveat sections evolve.

Add focused tests for the default report renderer:

```text
empty/no-region record
single region
multiple products
diagnostics
ClusterBlast present/absent
provenance rendering
special characters / Markdown escaping where applicable
deterministic repeated render
```

Use semantic substring assertions for flexible prose and exact golden snapshots only for intentionally stable machine-like blocks.

Avoid one giant brittle full-report snapshot unless the report is explicitly versioned as a stable format.

---

# 12. P8 — Lower priority: distinguish positional NRPS starter from unusual starter chemistry

This is a carry-over design limitation, not a dangerous bug.

Current mass eligibility blocks every:

```python
module.starter is True
```

because starter chemistry *may* include unresolved acyl/tail chemistry.

Upstream antiSMASH starter semantics are broader: an ordinary first loader-only NRPS module can be marked as a starter simply because it lacks a preceding C/KS domain.

Therefore the current gate can suppress otherwise valid peptide-backbone mass candidates.

That failure mode is conservative:

```text
valid candidate → null
```

rather than:

```text
wrong mass → emitted
```

so it should remain below provenance/schema work.

## Future audit

Classify starter evidence into something like:

```text
positional_starter
explicit_noncanonical_starter
unknown_starter
```

using audited domain evidence.

Potential explicit unusual starter signals include source-audited cases such as:

```text
Condensation_Starter
CAL_domain
SAT
other known nonstandard starter domains
```

Do not simply remove the starter gate globally.

Add real fixtures before broadening eligibility.

---

# 13. Suggested implementation order

```text
P0  antiSMASH 7.0.x result schema-2 support
    ├── parser allow-list
    ├── v7.0 fixture
    ├── rename v7.1 fixture
    └── correct support docs

P1  ClusterBlast data_version provenance
    ├── model field
    ├── parser validation
    ├── record schema 0.4.0
    └── exporter/tests

P2  fixture taxonomy and provenance
    ├── minimal parser-contract fixtures
    ├── serializer-faithful golden fixtures
    ├── README source matrix
    └── compatibility-only schema1 labeling

P3  diagnostic ownership
    ├── route duplicate diagnostics to actual record
    ├── add LoadedReviewInput.input_diagnostics
    ├── keep unknown parse failures input-scoped
    └── prevent record-comparison contamination

P4  isolated wheel/package CI smoke

P5  private/optional real-data integration regression harness

P6  dedicated Markdown renderer tests

P7  source-audit starter-module mass policy

THEN
    KnownClusterBlast coverage metrics
    LC-HRMS matching
    CI fail-on-new-diagnostic gates
    agent-context / analysis-digest export
```

P0–P3 are the most coherent batch for a `0.4.0` release because they all concern version-aware ClusterBlast compatibility and provenance.

---

# 14. Proposed `0.4.0` public contract

If P0–P3 are implemented together, a clean release boundary would be:

```text
package version                  0.4.0
record schema                    0.4.0
assemblyline schema              0.3.0
comparison schema                0.2.0 or 0.3.0*
cohort schema                    0.1.0
provenance manifest schema       0.1.0 or 0.2.0*
```

`*` Bump only if the corresponding serialized public shape actually changes.

Suggested record-schema migration note:

> Record schema 0.4.0 adds ClusterBlast `data_version` provenance and clarifies diagnostic ownership. Record-scoped diagnostics remain nested with their owning record; input-level enrichment diagnostics are represented separately rather than being assigned to the first parsed record.

No assembly-line schema bump is needed unless the assembly-line output shape or interpretation changes again.

---

# 15. Regression checklist

Before merging the follow-up patch, run:

```bash
ruff check antismash_review tests
ruff format --check antismash_review tests
mypy antismash_review
pytest --cov=antismash_review --cov-report=term-missing --cov-fail-under=90
python -m build
```

Required behavioral cases:

### ClusterBlast versions

```text
6.1 / result schema 1         PASS
7.0 / result schema 2         PASS
7.1 / result schema 3         PASS
8.0 / result schema 5         PASS
unsupported schema 4/6        FAIL
```

### Version provenance

```text
data_version present          preserved exactly
data_version absent           None
invalid data_version type     FAIL
section record_id mismatch    FAIL
section search_type mismatch  FAIL
```

### Diagnostic ownership

```text
parse failure, unknown record
    → input diagnostic, record_id=None

duplicate result for record B
    → B diagnostic only

unattachable result for known B
    → B diagnostic

unattachable result for unknown record
    → input diagnostic, retains source record_id if known

valid A + bad B + valid C in lenient mode
    → A and C remain attached
```

### Comparison

```text
input-level parse warning
    → must not appear as a new diagnostic for arbitrary first record

record-B duplicate warning
    → must not appear as a diagnostic delta for record A
```

### Packaging

```text
wheel built
fresh venv install
pip check passes
CLI works from /tmp
import path points to venv site-packages
```

---

# 16. What is already closed and should not be reopened casually

The following areas looked substantially improved at reviewed HEAD and should remain stable unless new evidence appears:

- free amino-acid chemistry and peptide dehydration;
- one-residue cyclic mass guard;
- module-level incorporation versus raw monomer-pairing multiplicity;
- conflicting versus identical duplicate monomer calls;
- cross-CDS duplicate integrity flags;
- incomplete/missing-domain mass gating;
- architecture scoring restricted to `nrps_pks_domains`;
- explicit NRPS module typing;
- evidence retention for unsupported product classes;
- transactional ClusterBlast attachment in strict/lenient modes;
- correct v6/v7.1/v8 outer module-container schema 2 in version-labeled fixtures;
- record JSON provenance object already embedded in schema 0.3.0;
- 90% CI coverage floor.

The next patch should avoid refactoring these stable areas unless required by P0–P3.

---

# 17. Final recommendation

The codebase is in a markedly better state than before the v0.3 series. The next step should **not** be a large feature expansion.

A focused `0.4.0` provenance/compatibility release would give the strongest return:

1. close the real antiSMASH 7.0 schema-2 hole;
2. preserve ClusterBlast `data_version`;
3. make fixture claims auditable and historically precise;
4. separate input-level diagnostics from record-level diagnostics;
5. then add isolated package and private real-output regression tests.

After those are complete, KnownClusterBlast-derived coverage metrics and LC-HRMS matching will sit on a much more defensible provenance foundation.

The key design principle should remain:

> **Preserve what antiSMASH actually emitted, record which antiSMASH/data generation emitted it, and keep uncertainty at the same scope as the evidence.**
