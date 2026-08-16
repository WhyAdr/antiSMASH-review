# antiSMASH-review Monomer-Pairing Mass Gate Patch Plan (`monomer-patch-plan.md`)

**Repository:** `WhyAdr/antiSMASH-review`\
**Target branch:** `main`\
**Primary scope:** `antismash_review/assemblyline.py`\
**Supporting scope:** models/exports/tests/schema/semantic contract\
**Motivating antiSMASH version:** 8.0.4\
**Proposed package/schema target:** `0.3.0`

---

## 1. Executive summary

The current assembly-line implementation correctly preserves antiSMASH `/monomer_pairings` as raw evidence, but it makes one important semantic mistake downstream:

> every `/monomer_pairings` entry is currently flattened into `AssemblyLinePrediction.chain` as though each entry represented an independent biosynthetic incorporation.

That assumption is not generally safe.

A real antiSMASH 8.0.4 example from `SM-ZPG19` demonstrates the failure mode:

```text
one aSModule
├── locus_tags
│   ├── Z1919F_012860
│   └── Z1919F_012855
├── domains
│   ├── Condensation_Starter
│   ├── AMP-binding
│   ├── PCP
│   ├── Epimerization
│   └── TIGR01720
└── monomer_pairings
    ├── Orn -> D-Orn
    └── Orn -> D-Orn
```

This is **one antiSMASH module**, not two consecutive D-Orn modules.

antiSMASH 8.0.4 appears to associate the same cross-CDS `Module` object with both constituent CDSs and subsequently append the same monomer prediction twice during NRPS/PKS result annotation. The GenBank exporter then faithfully serializes both stored values.

Consequently:

```text
raw annotation evidence:
    Orn -> D-Orn
    Orn -> D-Orn

must NOT automatically become:

biological chain:
    D-Orn -> D-Orn
```

Instead:

```text
raw annotation evidence:
    Orn -> D-Orn
    Orn -> D-Orn

interpreted incorporation:
    D-Orn

structural signal:
    cross-CDS duplicate monomer annotation
```

The patch should therefore introduce a strict distinction between:

1. **raw monomer evidence**, and
2. **module-level interpreted incorporation**.

Raw evidence must remain lossless.

Mass estimation and biological chain construction must operate on interpreted module-level incorporation slots rather than raw pairing multiplicity.

The same behavior can also be exploited as a useful structural-integrity review signal. A duplicated pairing inside a multi-CDS `aSModule` can identify loci worth inspecting for gene fragmentation, pseudogenization, gene-calling artifacts, or assembly errors.

It must **not**, however, be reported as proof of any of those states.

---

# 2. Current behavior

The parser currently does the correct thing at ingestion time.

`antismash_review/genbank.py::_module()` preserves:

```python
monomer_pairings=_values(q, "monomer_pairings")
```

without deduplicating or interpreting them.

This behavior should remain unchanged.

The problem enters later in `antismash_review/assemblyline.py`.

Currently:

```python
if module.monomer_pairings:
    calls = tuple(_pairing_call(raw) for raw in module.monomer_pairings)
```

and later:

```python
chain = tuple(
    call
    for module in module_predictions
    for call in module.monomer_calls
)
```

Thus:

```text
Module 1:
    Orn -> D-Orn
    Orn -> D-Orn
```

becomes:

```text
chain = (D-Orn, D-Orn)
```

despite there being only one `aSModule`.

`_estimate_core_mass()` then consumes `chain`, so annotation multiplicity can become monomer stoichiometry.

This violates the intended semantic hierarchy:

```text
aSModule
    ↓
biosynthetic module identity
    ↓
one default collinear incorporation slot

/monomer_pairings
    ↓
evidence attached to that module
```

The parser should never silently replace the first relationship with the second.

---

# 3. Core semantic rule

Adopt the following invariant:

> **Module identity determines the default incorporation slot. Pairing multiplicity describes evidence associated with that module.**

For non-iterative collinear interpretation:

```text
1 aSModule ≈ at most 1 interpreted incorporation slot
```

This does **not** mean an NRPS module can never be reused biochemically.

Iterative/non-collinear NRPS behavior exists.

It means only that:

> multiple strings in `/monomer_pairings` do not independently prove repeated catalytic turnover.

Iteration must be represented as a separate uncertainty about **incorporation count**.

---

# 4. Preserve two semantic layers

The patch should explicitly distinguish:

## Layer A — raw monomer evidence

Keep:

```python
ModulePrediction.monomer_calls
```

as the complete parsed representation of all antiSMASH-derived calls.

Example:

```python
monomer_calls = (
    MonomerCall(
        substrate="Orn",
        monomer="D-Orn",
        ...
    ),
    MonomerCall(
        substrate="Orn",
        monomer="D-Orn",
        ...
    ),
)
```

Nothing is discarded.

---

## Layer B — interpreted module incorporation

Add one derived module-level call representing the default collinear biological interpretation.

Suggested field:

```python
ModulePrediction.incorporation_call
```

For the SM-ZPG19 case:

```python
incorporation_call = MonomerCall(
    substrate="Orn",
    monomer="D-Orn",
    display="D-Orn",
    source="module_pairing",
    confidence="high",
    notes=(
        "2 identical raw module pairings collapsed to one incorporation slot",
    ),
)
```

Then:

```python
AssemblyLinePrediction.chain
```

must be built from:

```python
module.incorporation_call
```

rather than:

```python
module.monomer_calls
```

---

# 5. Proposed `ModulePrediction` additions

Extend the dataclass approximately as follows:

```python
PairingStatus = Literal[
    "single",
    "identical_duplicate",
    "conflicting",
    "specificity_fallback",
    "unresolved",
]
```

and:

```python
@dataclass(slots=True, frozen=True)
class ModulePrediction:
    ...
    monomer_calls: tuple[MonomerCall, ...]   # raw evidence
    incorporation_call: MonomerCall          # interpreted slot

    pairing_status: PairingStatus
    raw_pairing_count: int
    unique_pairing_count: int

    integrity_flags: tuple[str, ...]
    ...
```

`monomer_calls` should deliberately retain its existing meaning as raw evidence to minimize API churn.

`incorporation_call` becomes the value used for:

- `AssemblyLinePrediction.chain`
- resolved monomer counts
- peptide-core stoichiometry
- mass estimation

---

# 6. Pairing interpretation rules

## Case A — one pairing

Input:

```text
Ser -> Ser
```

Result:

```text
raw calls:        [Ser]
incorporation:    Ser
pairing_status:   single
```

No behavior change.

---

## Case B — identical duplicate calls in one module

Input:

```text
Orn -> D-Orn
Orn -> D-Orn
```

Result:

```text
raw calls:
    D-Orn
    D-Orn

incorporation:
    D-Orn

pairing_status:
    identical_duplicate

raw_pairing_count:
    2

unique_pairing_count:
    1
```

The raw duplicate remains visible.

The biological chain receives only one D-Orn slot.

---

## Case C — identical duplicates in a multi-CDS module

Input:

```text
aSModule
locus_tags:
    CDS_A
    CDS_B

monomer_pairings:
    Orn -> D-Orn
    Orn -> D-Orn
```

Result:

```text
incorporation:
    D-Orn

integrity_flags:
    duplicate_monomer_pairing
    cross_cds_duplicate_monomer_pairing
```

Also emit an evidence-scoped warning/caveat such as:

> Identical monomer pairing evidence occurs multiple times within one multi-CDS antiSMASH module. One incorporation slot is used for assembly-line interpretation; raw duplicates are retained. The cross-CDS pattern may warrant structural-integrity review but does not prove pseudogenization or assembly error.

---

## Case D — multiple different pairings in one module

Input:

```text
Ser -> Ser
Leu -> Leu
```

This must **not** become:

```text
Ser -> Leu
```

because there is still only one module.

Instead:

```text
raw calls:
    Ser
    Leu

incorporation:
    ?

pairing_status:
    conflicting

confidence:
    unresolved
```

The `incorporation_call.notes` should preserve the alternatives:

```text
conflicting module pairing evidence:
Ser -> Ser | Leu -> Leu
```

`ModulePrediction.monomer_calls` retains both typed alternatives.

Mass estimation remains unavailable.

---

## Case E — two separate modules with identical monomers

Input:

```text
Module 1:
    Orn -> D-Orn

Module 2:
    Orn -> D-Orn
```

Result:

```text
chain:
    D-Orn -> D-Orn
```

This is a genuine two-slot assembly-line interpretation.

**Never deduplicate across module boundaries.**

The scope of duplicate collapse is strictly:

```text
within the same aSModule
```

not:

```text
within the predicted chain
```

---

# 7. Comparison key for duplicate calls

Do not use a global `set()` on raw strings.

Raw strings must remain untouched.

For classification, compare parsed calls using a conservative key such as:

```python
def _call_identity(call: MonomerCall) -> tuple[str | None, str | None]:
    return call.substrate, call.monomer
```

Because `_pairing_call()` already strips surrounding whitespace, this should be sufficient initially.

Do not aggressively normalize chemical tokens.

For example, avoid transformations that could accidentally erase meaningful distinctions involving:

- stereochemistry,
- methylation,
- oxidation state,
- unusual extender-unit syntax.

Preservation is more important than cosmetic canonicalization.

---

# 8. Suggested helper architecture

Introduce a dedicated interpreter rather than embedding everything inside `_calls_for_module()`.

For example:

```python
def _interpret_module_calls(
    module: Module,
    calls: tuple[MonomerCall, ...],
) -> tuple[
    MonomerCall,
    PairingStatus,
    tuple[str, ...],
]:
    ...
```

or preferably return a small typed object:

```python
@dataclass(slots=True, frozen=True)
class ModuleCallInterpretation:
    incorporation_call: MonomerCall
    pairing_status: PairingStatus
    raw_pairing_count: int
    unique_pairing_count: int
    integrity_flags: tuple[str, ...]
    warnings: tuple[str, ...]
```

This isolates three concerns:

```text
parsing
    ↓
raw MonomerCall objects

interpretation
    ↓
one module incorporation slot

assembly-line construction
    ↓
ordered incorporation chain
```

That separation will make future handling of non-collinear NRPS behavior much safer.

---

# 9. `_calls_for_module()` responsibility

Keep `_calls_for_module()` focused primarily on **evidence retrieval**:

```text
module pairing
    >
domain specificity fallback
    >
unknown
```

It should not be responsible for chain stoichiometry.

Recommended flow:

```python
raw_calls, warnings = _calls_for_module(module, domains)

interpretation = _interpret_module_calls(
    module,
    raw_calls,
)

return ModulePrediction(
    ...
    monomer_calls=raw_calls,
    incorporation_call=interpretation.incorporation_call,
    pairing_status=interpretation.pairing_status,
    raw_pairing_count=interpretation.raw_pairing_count,
    unique_pairing_count=interpretation.unique_pairing_count,
    integrity_flags=interpretation.integrity_flags,
    warnings=(
        *warnings,
        *interpretation.warnings,
    ),
)
```

---

# 10. Fix `predict_assembly_lines()`

Replace:

```python
chain = tuple(
    call
    for module in module_predictions
    for call in module.monomer_calls
)
```

with:

```python
chain = tuple(
    module.incorporation_call
    for module in module_predictions
)
```

assuming every module receives one interpreted slot.

This makes `chain` explicitly:

> an ordered sequence of module-level incorporation hypotheses

rather than:

> a flattened list of annotation strings.

Update the docstring accordingly.

---

# 11. Mass-estimator correction

`_estimate_core_mass()` can remain structurally similar once `chain` is corrected.

Its semantics should become:

```text
total_monomers
    =
number of interpreted incorporation slots
```

not:

```text
number of raw monomer_pairing qualifiers
```

For SM-ZPG19-like evidence:

```text
raw calls:              2
unique raw calls:       1
incorporation slots:    1
```

Therefore:

```text
MassEstimate.total_monomers == 1
```

rather than `2`.

---

# 12. Conflicting calls and mass coverage

For:

```text
one module
two conflicting pairing calls
```

the interpreted chain should contain:

```text
?
```

once.

Therefore:

```text
total_monomers = 1
resolved_monomers = 0
coverage_fraction = 0.0
```

rather than:

```text
total_monomers = 2
```

This more accurately represents the biological uncertainty:

> one incorporation site exists, but its substrate identity is unresolved.

---

# 13. Iterative-module semantics

Do not treat repeated `/monomer_pairings` as evidence of iteration.

Keep antiSMASH's existing `/iterative` flag as upstream evidence, but interpret it conservatively.

The semantic contract should explicitly state:

> The antiSMASH `iterative` flag is retained as upstream module evidence. `antiSMASH-review` does not infer the number of catalytic turnovers from that flag and does not treat absence of the flag as proof that an NRPS module cannot behave iteratively.

For an iterative module:

```text
known monomer:
    Ser

known module:
    iterative=True
```

the parser may still represent:

```text
incorporation_call = Ser
```

as the module's substrate identity.

But:

```text
incorporation count = unresolved
```

and core-mass estimation must remain blocked.

The chain is therefore a **module-level sequence hypothesis**, not necessarily complete final-product stoichiometry for iterative systems.

Do not expand:

```text
Ser
```

into:

```text
Ser -> Ser -> Ser
```

without independent evidence for turnover count.

---

# 14. Cross-CDS duplicate calls as an integrity signal

This antiSMASH behavior should not merely be “fixed and forgotten.”

Preserve it as a potentially useful QC signal.

Recommended flags:

```text
duplicate_monomer_pairing
cross_cds_duplicate_monomer_pairing
```

The second is particularly interesting because it can prioritize loci for manual inspection.

Potential biological explanations include:

1. legitimate naturally split NRPS/PKS architecture;
2. gene-calling fragmentation;
3. true disruptive mutation;
4. pseudogenization;
5. assembly or polishing error;
6. frameshift;
7. premature stop;
8. indel within an otherwise conserved biosynthetic gene.

The flag itself must remain neutral.

Recommended human-facing wording:

> Cross-CDS duplicate monomer evidence detected. This module may merit inspection for gene fragmentation, pseudogenization, or assembly/gene-calling artifacts; the annotation pattern alone does not establish any of these states.

---

# 15. Do not automatically emit `pseudogene=true`

This patch must not make:

```text
cross-CDS module
+
duplicate pairing
```

equivalent to:

```text
pseudogene
```

The existing semantic contract correctly distinguishes evidence from biological conclusion.

Maintain that philosophy.

A future high-confidence pseudogenization assessment should require independent evidence such as:

```text
explicit /pseudo or /pseudogene annotation
        +
sequence disruption
        +
comparative homolog evidence
```

or equivalent multi-layer support.

---

# 16. Future structural-integrity classifier

This patch can establish the evidence hooks without implementing a full pseudogene classifier.

A later module could distinguish:

```text
cross_cds_module
        ↓
structural integrity assessment
        ├── likely bona-fide split
        ├── possible gene-call artifact
        ├── possible assembly artifact
        ├── fragmentation candidate
        └── likely pseudogenized
```

Potential evidence layers:

### Layer 1 — antiSMASH structure

- multi-CDS `aSModule`
- duplicate pairing within module
- domain distribution across parent CDSs
- incomplete/complete flag
- unexpected domain boundary

### Layer 2 — GenBank annotation

- `/pseudo`
- `/pseudogene`
- partial/fuzzy CDS
- internal stop-related notes
- adjacent CDS arrangement
- short intergenic separation

### Layer 3 — comparative genomics

Compare close homologs:

```text
reference strains:
    C-A-PCP-E = one protein

query strain:
    C-A | A-PCP-E = two proteins
```

Strong conservation of an intact single-CDS homolog would substantially strengthen fragmentation suspicion.

### Layer 4 — nucleotide evidence

Inspect:

- stop codons,
- frameshift indels,
- interrupted reading frame,
- missing conserved segment,
- splice-like or translational anomalies.

### Layer 5 — raw sequencing evidence

For ONT or other long-read assemblies:

- remap raw reads;
- inspect indel support;
- examine homopolymer contexts;
- compare polishing stages.

This is the level required to distinguish:

```text
real pseudogenization
```

from:

```text
one-base polishing error
```

with meaningful confidence.

---

# 17. SM-ZPG19 motivating example

Document the observed pattern as a regression case.

Conceptually:

```text
SM-ZPG19

one aSModule
┌─────────────────────────────────────────────┐
│ CDS Z1919F_012860                           │
│     Condensation_Starter                    │
│                   ↓                         │
│ CDS Z1919F_012855                           │
│     AMP-binding → PCP → Epimerization       │
└─────────────────────────────────────────────┘

raw antiSMASH calls:
    Orn -> D-Orn
    Orn -> D-Orn

current antiSMASH-review interpretation:
    D-Orn -> D-Orn        [wrong]

patched interpretation:
    D-Orn                 [module-level slot]

additional signal:
    cross_cds_duplicate_monomer_pairing
```

The intact PA-LS101 homolog provides a useful conceptual control:

```text
one CDS
Cstarter → A(Orn) → PCP → E
                    ↓
                  D-Orn
```

This comparative observation motivated the patch but should **not** be hard-coded into generic parser logic.

---

# 18. Export changes

## JSON

JSON already serializes dataclasses via `asdict()`, so new fields should naturally appear.

Example:

```json
{
  "monomer_calls": [
    {
      "substrate": "Orn",
      "monomer": "D-Orn"
    },
    {
      "substrate": "Orn",
      "monomer": "D-Orn"
    }
  ],
  "incorporation_call": {
    "substrate": "Orn",
    "monomer": "D-Orn"
  },
  "pairing_status": "identical_duplicate",
  "raw_pairing_count": 2,
  "unique_pairing_count": 1,
  "integrity_flags": [
    "duplicate_monomer_pairing",
    "cross_cds_duplicate_monomer_pairing"
  ]
}
```

This is the preferred machine-readable representation.

---

## Markdown

Change the module table from approximately:

```text
| Module | Locus tags | Type | Calls | Domains | ... |
```

to:

```text
| Module | Locus tags | Type | Raw calls | Incorporation | Pairing status | Integrity flags | Domains | ... |
```

SM-ZPG19-like output:

```text
Raw calls:
D-Orn, D-Orn

Incorporation:
D-Orn

Pairing status:
identical_duplicate

Integrity:
cross_cds_duplicate_monomer_pairing
```

This makes the anomaly immediately visible without corrupting the biological chain.

---

## TSV

Preserve the current evidence-oriented behavior of one row per raw call if possible.

Do **not** silently convert TSV into one-row-per-module without a deliberate migration.

Add columns such as:

```text
raw_call_index
raw_pairing_count
unique_pairing_count
pairing_status
interpreted_substrate
interpreted_monomer
interpreted_call_confidence
integrity_flags
```

Thus two duplicate raw rows remain visible:

```text
raw row 1: D-Orn
raw row 2: D-Orn
```

while both rows clearly report:

```text
interpreted_monomer = D-Orn
raw_pairing_count = 2
unique_pairing_count = 1
pairing_status = identical_duplicate
```

This preserves forensic transparency.

---

# 19. Assembly-line schema bump

The current assembly-line schema is `0.2.0`.

This patch changes the semantics of:

```text
AssemblyLinePrediction.chain
MassEstimate.total_monomers
MassEstimate.resolved_monomers
coverage_fraction
```

for modules containing multiple pairing calls.

It also adds new module-level fields.

Therefore bump:

```python
ASSEMBLYLINE_SCHEMA_VERSION = "0.3.0"
```

Recommended package bump:

```python
__version__ = "0.3.0"
```

Do not change the record schema merely because derived assembly-line output changes.

The parsed `Record` representation remains lossless and unchanged.

---

# 20. Regeneration notice

Add a migration note to `references/semantic-contract.md`.

Suggested wording:

> **Assembly-line monomer multiplicity notice:** assembly-line outputs produced by parser version `0.2.0` interpreted every `/monomer_pairings` value as a separate chain entry. antiSMASH outputs containing multiple pairings within one `aSModule`, particularly duplicated calls on cross-CDS modules, could therefore overcount predicted incorporation slots and distort derived core-mass candidates. Regenerate affected assembly-line outputs with `parser_version >= 0.3.0`. Raw parsed `Module.monomer_pairings` evidence itself was preserved correctly and is not affected.

This should sit alongside the existing `0.1.0` core-mass regeneration notice.

---

# 21. Semantic-contract updates

Update the assembly-line section to state explicitly:

### Raw evidence

> `Module.monomer_pairings` are retained in source order without deduplication.

### Interpretation

> Pairing multiplicity within one `aSModule` does not by itself imply multiple incorporated residues.

### Duplicate calls

> Identical repeated pairings within one module are retained as raw evidence but collapse to one default collinear incorporation slot.

### Conflicting calls

> Distinct pairings within one module represent conflicting/alternative evidence and produce one unresolved incorporation slot rather than multiple consecutive residues.

### Module boundaries

> Identical monomers belonging to distinct `aSModule` objects remain distinct incorporation slots.

### Cross-CDS signal

> Duplicate pairings on a multi-CDS module may be exposed as a structural-integrity review signal but do not prove fragmentation, pseudogenization, frameshift, or assembly error.

### Iteration

> The upstream `iterative` flag is preserved but does not specify repeat count and is not treated as a general proof or disproof of iterative NRPS chemistry.

---

# 22. Tests

Add focused tests to `tests/test_assemblyline.py`.

## Test 1 — identical duplicate within one module

```python
pairings=[
    "Orn -> D-Orn",
    "Orn -> D-Orn",
]
```

Assert:

```python
len(module.monomer_calls) == 2

module.incorporation_call.monomer == "D-Orn"

prediction.chain == (D-Orn,)

module.pairing_status == "identical_duplicate"

module.raw_pairing_count == 2
module.unique_pairing_count == 1
```

---

## Test 2 — cross-CDS duplicate

Use:

```python
locus_tags=["CDS_A", "CDS_B"]
```

Assert:

```python
module.multi_cds is True

"cross_cds_duplicate_monomer_pairing"
    in module.integrity_flags
```

Assert chain length is `1`.

---

## Test 3 — duplicate affects mass count correctly

Use a non-starter synthetic NRPS module with a modeled amino acid.

Before patch:

```text
total_monomers = 2
```

After patch:

```text
total_monomers = 1
```

Assert no duplicated dehydration event occurs.

---

## Test 4 — conflicting pairings

Input:

```python
pairings=[
    "Ser -> Ser",
    "Leu -> Leu",
]
```

Assert:

```python
len(module.monomer_calls) == 2

len(prediction.chain) == 1

prediction.chain[0].monomer is None
prediction.chain[0].confidence == "unresolved"

module.pairing_status == "conflicting"

mass.linear_core_mass_da is None
```

---

## Test 5 — two modules with same monomer

Input:

```text
Module 1:
    Orn -> D-Orn

Module 2:
    Orn -> D-Orn
```

Assert:

```python
prediction.chain == (
    D-Orn,
    D-Orn,
)
```

This prevents an overly broad deduplication regression.

---

## Test 6 — iterative module

Assert:

```python
module.iterative is True
```

does not automatically expand one pairing into repeated residues.

Mass must remain unavailable because incorporation count is unresolved.

---

## Test 7 — specificity fallback

Ensure existing behavior survives:

```text
no /monomer_pairings
AMP-binding specificity:
    substrate consensus: Ser
```

Result:

```text
one low-confidence incorporation slot
```

No regression.

---

## Test 8 — malformed pairing

A malformed raw call remains present in:

```python
module.monomer_calls
```

while:

```python
incorporation_call
```

remains one unresolved slot.

---

## Test 9 — Markdown export

Assert output exposes:

```text
Raw calls
Incorporation
Pairing status
Integrity flags
```

---

## Test 10 — TSV export

Verify duplicate raw rows remain present while interpreted monomer is identical and module-level count remains one.

---

## Test 11 — JSON schema

Assert:

```json
"schema_version": "0.3.0"
```

and presence of:

```text
incorporation_call
pairing_status
raw_pairing_count
unique_pairing_count
integrity_flags
```

---

# 23. Integration fixture

Prefer a compact synthetic antiSMASH-like GenBank regression fixture over committing an entire real bacterial region.

Extend `tests/fixtures/build_fixture.py` or add a dedicated helper producing:

```text
CDS_A
    Condensation_Starter

CDS_B
    AMP-binding
    PCP
    Epimerization

aSModule
    locus_tags:
        CDS_A
        CDS_B

    monomer_pairings:
        Orn -> D-Orn
        Orn -> D-Orn
```

The fixture should demonstrate only the representation needed to reproduce the bug.

Optionally construct a paired control:

```text
CDS_INTACT
    Cstarter
    AMP-binding
    PCP
    Epimerization

monomer_pairings:
    Orn -> D-Orn
```

Expected interpreted output from both:

```text
D-Orn
```

with the split fixture additionally carrying:

```text
cross_cds_duplicate_monomer_pairing
```

Real SM-ZPG19 and PA-LS101 files can remain local/manual validation cases unless repository size or provenance requirements justify committing them.

---

# 24. Suggested implementation sequence

## Patch 1 — representation

Modify:

```text
antismash_review/assemblyline.py
```

Add:

- `PairingStatus`
- module-call interpretation helper
- `incorporation_call`
- pairing counts
- integrity flags

Do not modify GenBank parsing.

---

## Patch 2 — chain semantics

Change:

```python
chain = flatten(module.monomer_calls)
```

to:

```python
chain = one incorporation_call per module
```

Run focused assembly-line tests.

---

## Patch 3 — mass regression

Verify `_estimate_core_mass()` consumes only interpreted chain slots.

Add duplicate and conflicting-call mass tests.

---

## Patch 4 — exporters

Update:

```text
antismash_review/exporters/assemblyline_markdown.py
antismash_review/exporters/assemblyline_table.py
```

JSON should largely follow automatically through `asdict()`.

---

## Patch 5 — schema/version

Update:

```text
antismash_review/schema.py
antismash_review/_version.py
```

Target:

```text
0.3.0
```

---

## Patch 6 — semantic contract

Update:

```text
references/semantic-contract.md
```

Document:

- raw/effective distinction;
- module-scoped deduplication;
- conflicting calls;
- iterative uncertainty;
- cross-CDS review signal;
- regeneration notice.

---

## Patch 7 — full regression

Run:

```bash
pytest -q
```

Then preferably:

```bash
pytest --cov=antismash_review --cov-report=term-missing
```

and verify deterministic exports remain byte-stable.

---

# 25. Acceptance criteria

The patch is complete when all of the following are true:

- [ ] GenBank parsing preserves duplicate `/monomer_pairings` exactly.
- [ ] `ModulePrediction.monomer_calls` retains every raw interpreted call.
- [ ] Each ordinary `aSModule` contributes at most one default incorporation slot.
- [ ] Identical duplicate calls inside one module collapse only at the interpretation layer.
- [ ] Distinct calls inside one module become one unresolved incorporation slot.
- [ ] Identical calls in different modules remain separate residues.
- [ ] `AssemblyLinePrediction.chain` represents module-level incorporation slots.
- [ ] `_estimate_core_mass()` uses the interpreted chain.
- [ ] Cross-CDS duplicate calls generate an explicit structural-integrity review signal.
- [ ] The signal does not claim pseudogenization or assembly error.
- [ ] Iterative modules do not receive an invented repeat count.
- [ ] JSON exposes both raw evidence and interpreted incorporation.
- [ ] Markdown clearly distinguishes raw calls from incorporation.
- [ ] TSV preserves raw duplicate evidence.
- [ ] Assembly-line schema is bumped to `0.3.0`.
- [ ] A regeneration notice documents the `0.2.0` multiplicity issue.
- [ ] Synthetic regression tests cover the SM-ZPG19-like edge case.
- [ ] Existing single-pairing behavior remains unchanged.
- [ ] Existing tests pass.

---

# 26. Expected behavior matrix

| Situation                             |                Raw calls |                                 Interpreted slots | Mass behavior     | Integrity signal    |
| ------------------------------------- | -----------------------: | ------------------------------------------------: | ----------------- | ------------------- |
| One module, one `Ser` call            |                        1 |                                           1 × Ser | normal gate       | none                |
| One module, `D-Orn` duplicated        |                        2 |                                         1 × D-Orn | count as one      | duplicate pairing   |
| Multi-CDS module, `D-Orn` duplicated  |                        2 |                                         1 × D-Orn | count as one      | cross-CDS duplicate |
| One module, `Ser` + `Leu` conflict    |                        2 |                                    1 × unresolved | mass null         | conflicting pairing |
| Two modules, both `D-Orn`             |                  2 total |                                         2 × D-Orn | genuine diresidue | none                |
| Iterative module, one `Ser` call      |                        1 | one module-level Ser hypothesis; turnover unknown | mass null         | iterative caveat    |
| No pairing, one specificity consensus |               1 fallback |                             1 low-confidence slot | mass null         | fallback            |
| No usable evidence                    | 1 unresolved placeholder |                                 1 unresolved slot | mass null         | unresolved          |

---

# 27. Design principle

The most important architectural rule for this patch is:

```text
raw annotation multiplicity
        ≠
biosynthetic stoichiometry
```

Instead:

```text
antiSMASH aSModule
        ↓
module-level incorporation hypothesis
        ↓
assembly-line chain
        ↓
optional chemistry
```

while separately:

```text
raw qualifiers
        ↓
annotation anomalies
        ↓
cross-CDS / duplicate evidence
        ↓
structural-integrity review
```

This keeps the parser conservative while turning an antiSMASH edge case into useful biological QC information.

The result should be both **more chemically correct** and **more informative** than simply calling `set(monomer_pairings)` and discarding the duplication.
