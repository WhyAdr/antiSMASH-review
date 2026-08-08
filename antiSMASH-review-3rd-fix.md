# antiSMASH-review third-pass fix plan

This document records the post-refactor audit of commit `cc68c88` against
[`antiSMASH-review-2nd-refactor.md`](antiSMASH-review-2nd-refactor.md). The refactor is
substantial and mostly functional, but it is not ready for final sign-off: three behavioral
defects, a coverage shortfall, and documentation drift remain.

The snippets below are implementation-ready guidance. They have not been applied to the
package code.

## Audit outcome

| Area | Result | Evidence |
|---|---|---|
| Ruff | Pass | `python -m ruff check .` |
| Formatting | Pass | 27 files already formatted |
| Mypy | Pass | 17 source files checked |
| Tests | Pass | 81 passed |
| Coverage | **Fail** | 92%; the plan requires at least 93% |
| SM-ZPG19 private tests | Pass | 6 passed without skips |
| Aggregate ClusterBlast attachment | Pass | 1 record with exactly 30 text results |
| Region-only ClusterBlast attachment | Pass | 10 records with exactly 3 results each |
| CLI smoke tests | Pass | All 10 planned commands returned status 0 |
| Markdown artifact checks | Pass | UTF-8, fences, links, whitespace, and mojibake clean |
| Working tree after audit | Pass | Clean; private data remains ignored |
| Wheel/import/entry-point audit | Blocked | Host Python lacks `setuptools`, `wheel`, and `build` |

## Fix order

1. Close the input-overwrite hole.
2. Normalize every recognized ClusterBlast JSON failure to `ClusterBlastParseError`.
3. Stop using record IDs as internal identities in coordinate matching.
4. Add the missing regression matrix and restore coverage to at least 93%.
5. Synchronize output wording and documentation.
6. Repeat the complete release gate in a build-capable environment.

## 1. Protect every discovered GenBank input

**Severity: high**

### Finding

[`antismash_review/cli.py`](antismash_review/cli.py) constructs `all_inputs` from `paths`,
which contains either aggregate GenBanks or region GenBanks. If both representations are
present, the unselected representation is omitted from overwrite protection.

A controlled reproduction used one aggregate GBK and one region GBK in the same temporary
directory, then selected the region GBK as `--output`. The command returned status 0 and
replaced the region GBK with Markdown beginning `# antiSMASH review`.

This directly violates the second-refactor requirement that overwrite protection cover every
GenBank and JSON input, not only files selected for parsing.

### Implementation diff

```diff
diff --git a/antismash_review/cli.py b/antismash_review/cli.py
@@
 def _load_review_records(
     manifest: InputManifest,
     *,
     lenient: bool,
 ) -> tuple[list[Record], set[Path]]:
     paths = manifest.aggregate_genbanks or manifest.region_genbanks
     all_inputs = (
-        set(paths)
+        set(manifest.aggregate_genbanks)
+        | set(manifest.region_genbanks)
         | set(manifest.json_files)
         | set(manifest.clusterblast_files)
         | set(manifest.knownclusterblast_files)
         | set(manifest.subclusterblast_files)
     )
```

### Regression test diff

```diff
diff --git a/tests/test_cli.py b/tests/test_cli.py
@@
+def test_cli_refuses_to_overwrite_unselected_region_genbank(
+    tmp_path: Path,
+    capsys: object,
+) -> None:
+    aggregate = tmp_path / "sample.gbk"
+    region = tmp_path / "sample.region001.gbk"
+    shutil.copy(FIXTURE, aggregate)
+    shutil.copy(FIXTURE, region)
+    original = region.read_bytes()
+
+    assert main(["inspect", str(tmp_path), "--output", str(region)]) == 2
+    captured = capsys.readouterr()  # type: ignore[union-attr]
+    assert "refusing to overwrite input file" in captured.err
+    assert region.read_bytes() == original
```

Apply the same assertion through `compare` with the vulnerable path appearing on either side.
The shared loader should make the code fix common to both commands.

### Acceptance

- Aggregate preference remains unchanged.
- Output cannot overwrite an aggregate GBK, a region GBK, or native JSON discovered on either
  comparison side.
- Rejected output leaves the input byte-for-byte unchanged and returns status 2.

## 2. Make ClusterBlast JSON failures public and strict

**Severity: high**

### Finding

[`antismash_review/clusterblast.py`](antismash_review/clusterblast.py) currently has two
failure-path violations:

- A present ClusterBlast section with a non-dictionary value is silently skipped.
- Missing fields such as `region_number` can raise raw `KeyError`; other malformed shapes can
  leak `TypeError`, `ValueError`, or `AttributeError`.

Observed probes:

```text
nondict_section status=0 diagnostics=[...]  # no clusterblast_parse_failed
missing_region_number uncaught=KeyError: 'region_number'
```

The parser's public contract requires `ClusterBlastParseError`. In lenient CLI mode that error
must become an explicit `clusterblast_parse_failed` diagnostic; it must never be discarded or
escape as an implementation exception.

### Implementation diff

Rename the existing implementation to a private unchecked helper, validate structural
containers that are currently skipped, and keep one public exception boundary:

```diff
diff --git a/antismash_review/clusterblast.py b/antismash_review/clusterblast.py
@@
-def parse_clusterblast_json(path: Path) -> list[ClusterBlastResult]:
+def _parse_clusterblast_json_unchecked(path: Path) -> list[ClusterBlastResult]:
     """Parse only ClusterBlast modules from a native antiSMASH JSON sidecar."""
@@
-    for record_idx, rec in enumerate(data.get("records", [])):
+    raw_records = data["records"]
+    if not isinstance(raw_records, list):
+        raise ClusterBlastParseError(f"antiSMASH JSON records is not a list in {path}")
+
+    for record_idx, rec in enumerate(raw_records):
+        if not isinstance(rec, dict):
+            raise ClusterBlastParseError(
+                f"antiSMASH JSON record {record_idx} is not an object in {path}"
+            )
         rec_id = rec.get("id")
         modules = rec.get("modules", {})
+        if not isinstance(modules, dict):
+            raise ClusterBlastParseError(
+                f"antiSMASH JSON modules for record {record_idx} is not an object in {path}"
+            )
@@
             section = cb_mod[sec_key]
             if not isinstance(section, dict):
-                continue
+                raise ClusterBlastParseError(
+                    f"ClusterBlast {sec_key} section is not an object in {path}"
+                )
@@
     return results
+
+
+def parse_clusterblast_json(path: Path) -> list[ClusterBlastResult]:
+    """Parse supported ClusterBlast JSON or raise one public enrichment error."""
+    path = Path(path)
+    try:
+        return _parse_clusterblast_json_unchecked(path)
+    except ClusterBlastParseError:
+        raise
+    except (KeyError, TypeError, ValueError, AttributeError) as exc:
+        raise ClusterBlastParseError(
+            f"Malformed ClusterBlast JSON structure in {path}: {exc}"
+        ) from exc
```

Do not broaden this into native JSON-to-`Record` support. This remains a version-gated parser
for `records[*].modules["antismash.modules.clusterblast"]` only.

### Regression test diff

```diff
diff --git a/tests/test_clusterblast.py b/tests/test_clusterblast.py
@@
+@pytest.mark.parametrize(
+    "section",
+    [
+        [],
+        {"schema_version": 5, "results": [{}]},
+        {"schema_version": 5, "results": "not-a-list"},
+    ],
+)
+def test_parse_clusterblast_json_malformed_shapes_use_public_error(
+    tmp_path: Path,
+    section: object,
+) -> None:
+    path = tmp_path / "malformed.json"
+    path.write_text(
+        json.dumps(
+            {
+                "records": [
+                    {
+                        "id": "contig_1",
+                        "modules": {
+                            "antismash.modules.clusterblast": {
+                                "schema_version": 2,
+                                "record_id": "contig_1",
+                                "general": section,
+                            }
+                        },
+                    }
+                ]
+            }
+        ),
+        encoding="utf-8",
+    )
+
+    with pytest.raises(ClusterBlastParseError):
+        parse_clusterblast_json(path)
```

Extend the existing lenient CLI test with malformed JSON, not only corrupt text:

```diff
diff --git a/tests/test_clusterblast.py b/tests/test_clusterblast.py
@@
+def test_cli_lenient_clusterblast_json_failure_is_diagnostic(
+    tmp_path: Path,
+    capsys: object,
+) -> None:
+    shutil.copy(FIXTURES / "semantics.gb", tmp_path / "semantics.gb")
+    (tmp_path / "semantics.json").write_text(
+        json.dumps(
+            {
+                "records": [
+                    {
+                        "id": "SEMANTICS.1",
+                        "modules": {
+                            "antismash.modules.clusterblast": {
+                                "schema_version": 2,
+                                "record_id": "SEMANTICS.1",
+                                "general": {"schema_version": 5, "results": [{}]},
+                            }
+                        },
+                    }
+                ]
+            }
+        ),
+        encoding="utf-8",
+    )
+
+    assert main(["inspect", str(tmp_path)]) == 2
+    capsys.readouterr()  # type: ignore[union-attr]
+
+    assert main(["inspect", str(tmp_path), "--lenient", "--format", "json"]) == 0
+    captured = capsys.readouterr()  # type: ignore[union-attr]
+    document = json.loads(captured.out)
+    assert "clusterblast_parse_failed" in {
+        diagnostic["code"] for diagnostic in document["diagnostics"]
+    }
```

Also add the planned cases that remain absent: module missing, unsupported result schema,
record/module ID mismatch, malformed query string, multiple JSON records, all three search
types, and invalid pairing structures.

### Acceptance

- A valid JSON document without the ClusterBlast module still returns an empty list.
- Every recognized malformed ClusterBlast structure raises `ClusterBlastParseError` with its
  source path.
- Strict CLI mode returns status 2 without a traceback.
- Lenient mode returns status 0 and emits `clusterblast_parse_failed`.
- No `KeyError`, `TypeError`, `ValueError`, or `AttributeError` escapes the public parser.

## 3. Use record positions as coordinate-match identities

**Severity: medium**

### Finding

[`antismash_review/compare.py`](antismash_review/compare.py) stores coordinate candidates in a
dictionary keyed by `left.record_id` and stores matched right records by `right.record_id`.
Coordinate mode does not require unique IDs, and antiSMASH region records commonly repeat a
contig ID.

A controlled case with two left records named `X`, two distinct coordinate spans, and two
uniquely corresponding right records incorrectly raised:

```text
ValueError: Non-one-to-one coordinate match: right record B matched by both X and X
```

The records have distinct object positions and a valid one-to-one solution; only their display
IDs repeat.

### Implementation diff

Use list positions for internal identity and retain record IDs only in messages and exports:

```diff
diff --git a/antismash_review/compare.py b/antismash_review/compare.py
@@
-        candidates_for_left: dict[str, list[tuple[Record, CoordinateMatchEvidence]]] = {}
-        for left in left_records:
-            candidates: list[tuple[Record, CoordinateMatchEvidence]] = []
-            for right in right_records:
+        candidates_for_left: list[
+            tuple[Record, list[tuple[int, Record, CoordinateMatchEvidence]]]
+        ] = []
+        for left in left_records:
+            candidates: list[tuple[int, Record, CoordinateMatchEvidence]] = []
+            for right_index, right in enumerate(right_records):
                 evidence = _calculate_overlap_evidence(left, right)
                 if (
                     evidence.left_overlap_fraction >= min_reciprocal_overlap
                     and evidence.right_overlap_fraction >= min_reciprocal_overlap
                 ):
-                    candidates.append((right, evidence))
-            candidates_for_left[left.record_id] = candidates
+                    candidates.append((right_index, right, evidence))
+            candidates_for_left.append((left, candidates))

-        matched_right_records: dict[str, str] = {}
-        for left in left_records:
-            cand_list = candidates_for_left[left.record_id]
+        matched_right_records: dict[int, str] = {}
+        for left, cand_list in candidates_for_left:
             if len(cand_list) == 0:
                 unmatched_left.append(left.record_id)
@@
             else:
-                right, evidence = cand_list[0]
-                if right.record_id in matched_right_records:
-                    prev_left = matched_right_records[right.record_id]
+                right_index, right, evidence = cand_list[0]
+                if right_index in matched_right_records:
+                    prev_left = matched_right_records[right_index]
                     raise ValueError(
                         f"Non-one-to-one coordinate match: right record {right.record_id} matched "
                         f"by both {prev_left} and {left.record_id}"
                     )
-                matched_right_records[right.record_id] = left.record_id
+                matched_right_records[right_index] = left.record_id
@@
-        matched_right_set = set(matched_right_records.keys())
-        for right in right_records:
-            if right.record_id not in matched_right_set:
+        matched_right_set = set(matched_right_records)
+        for right_index, right in enumerate(right_records):
+            if right_index not in matched_right_set:
                 unmatched_right.append(right.record_id)
```

This preserves the deliberately conservative rule that more than one above-threshold candidate
is ambiguous. It only prevents repeated labels from corrupting record identity.

### Regression test diff

```diff
diff --git a/tests/test_compare.py b/tests/test_compare.py
@@
+def test_coordinate_overlap_allows_repeated_display_ids() -> None:
+    left_a = _minimal_record("REPEATED")
+    left_a.genes.append(_make_gene("LA", 0, 100))
+    left_b = _minimal_record("REPEATED")
+    left_b.genes.append(_make_gene("LB", 200, 300))
+
+    right_a = _minimal_record("RIGHT")
+    right_a.genes.append(_make_gene("RA", 0, 100))
+    right_b = _minimal_record("RIGHT")
+    right_b.genes.append(_make_gene("RB", 200, 300))
+
+    result = compare_records(
+        [left_a, left_b],
+        [right_a, right_b],
+        left_input=Path("left.gb"),
+        right_input=Path("right.gb"),
+        match_method="coordinate_overlap",
+        assume_shared_coordinate_system=True,
+    )
+
+    assert len(result.matched) == 2
+    assert result.unmatched_left == []
+    assert result.unmatched_right == []
+    evidence = [item.coordinate_evidence for item in result.matched]
+    assert all(item is not None for item in evidence)
+    assert [item.overlap_bp for item in evidence if item is not None] == [100, 100]
```

### Acceptance

- Repeated left or right record IDs do not overwrite internal candidate state.
- Truly ambiguous candidates and true multiple-left-to-one-right assignments still fail.
- Match evidence belongs to the actual left/right objects used in the comparison.
- Output remains deterministic even when display IDs repeat.

## 4. Complete the promised comparison test matrix

**Severity: medium**

The current `test_compare_coordinate_overlap_below_threshold_and_ambiguity` checks only the
below-threshold case. Add separate tests whose names and bodies correspond to the planned
contract.

```diff
diff --git a/tests/test_compare.py b/tests/test_compare.py
@@
+def test_coordinate_overlap_rejects_equal_candidate_tie() -> None:
+    left = _minimal_record("LEFT")
+    left.genes.append(_make_gene("L", 0, 100))
+    right_a = _minimal_record("RIGHT_A")
+    right_a.genes.append(_make_gene("RA", 0, 100))
+    right_b = _minimal_record("RIGHT_B")
+    right_b.genes.append(_make_gene("RB", 0, 100))
+
+    with pytest.raises(ValueError, match="Ambiguous coordinate match"):
+        compare_records(
+            [left],
+            [right_a, right_b],
+            left_input=Path("left.gb"),
+            right_input=Path("right.gb"),
+            match_method="coordinate_overlap",
+            assume_shared_coordinate_system=True,
+        )
+
+
+def test_coordinate_overlap_rejects_multiple_left_to_one_right() -> None:
+    left_a = _minimal_record("LEFT_A")
+    left_a.genes.append(_make_gene("LA", 0, 100))
+    left_b = _minimal_record("LEFT_B")
+    left_b.genes.append(_make_gene("LB", 0, 100))
+    right = _minimal_record("RIGHT")
+    right.genes.append(_make_gene("R", 0, 100))
+
+    with pytest.raises(ValueError, match="Non-one-to-one coordinate match"):
+        compare_records(
+            [left_a, left_b],
+            [right],
+            left_input=Path("left.gb"),
+            right_input=Path("right.gb"),
+            match_method="coordinate_overlap",
+            assume_shared_coordinate_system=True,
+        )
+
+
+@pytest.mark.parametrize("match_method", ["record_id", "record_region", "single_record"])
+def test_assumption_flag_rejected_for_every_non_coordinate_mode(
+    match_method: str,
+) -> None:
+    sem = FIXTURES / "semantics.gb"
+    with pytest.raises(SystemExit) as exc:
+        main(
+            [
+                "compare",
+                str(sem),
+                str(sem),
+                "--match-by",
+                match_method,
+                "--assume-shared-coordinate-system",
+            ]
+        )
+    assert exc.value.code == 2
+
+
+def test_intergenic_summary_handles_cross_origin_cds() -> None:
+    record = parse_genbank(FIXTURES / "cross-origin.gb")[0]
+    summary = intergenic_summary(record)
+    assert summary.circular_wrap_included
+    assert summary.gap_count == 1
+    assert summary.total_bp == 60
```

Add focused tests for repeated diagnostic fingerprints, sorted unmatched IDs, empty
`single_record` inputs, missing-number `record_region` inputs, and the right-side duplicate-ID
case in `record_id` mode.

Coverage must be at least the 93% baseline without excluding new modules or weakening the
measurement configuration.

## 5. Synchronize wording and semantic documentation

**Severity: low**

The user-facing term is **record matching mode**, not identity mode or the less precise match
method. The semantic contract also lacks several limitations explicitly required by the plan.

### Output wording diff

```diff
diff --git a/antismash_review/exporters/compare_markdown.py b/antismash_review/exporters/compare_markdown.py
@@
-        f"- Match method: `{result.match_method}`",
+        f"- Record matching mode: `{result.match_method}`",
```

### Semantic contract diff

```diff
diff --git a/references/semantic-contract.md b/references/semantic-contract.md
@@
 | `missing_nrps_pks_architecture` | WARNING | Region products contain NRPS/PKS terms but no aSTool=nrps_pks_domains domains were parsed | Does not prove the annotation is wrong; the relevant domains may be in a different record |
+| `clusterblast_parse_failed` | WARNING (lenient) | A recognized ClusterBlast text or JSON sidecar could not be parsed or attached safely | Does not repair or silently discard the malformed sidecar; strict mode fails instead |
@@
 ## ClusterBlast sidecar enrichment

 - Parse text sidecars from canonical `clusterblast/`, `knownclusterblast/`, and `subclusterblast/` directories using natural region sorting (`_cN.txt`).
 - Parse native antiSMASH JSON `antismash.modules.clusterblast` modules validating module schema `2` and result schema `5`.
 - Precedence is per `(record_id, region_number, search_type)` key: text files are preferred when present for a key, while JSON fills remaining keys.
+- Retain valid empty results as negative evidence and retain source path, SHA-256, source format, and supported schema versions as provenance.
 - In strict mode, malformed or unattached sidecars raise `ClusterBlastParseError`. In lenient mode, errors emit a `clusterblast_parse_failed` diagnostic.
+- Retain `misc_feature` entries only as raw GenBank evidence; do not infer that they are ClusterBlast results.
@@
- Supported matching modes:
+- Supported record matching modes:
@@
   - `coordinate_overlap`: Match records by reciprocal overlap of feature spans. Requires explicit `--assume-shared-coordinate-system` and a reciprocal overlap fraction in `(0, 1]` (default 0.80).
+- The coordinate assumption is appropriate only when coordinate correspondence has been independently established, such as re-annotations of the same contigs. It does not establish sequence homology between arbitrary isolates.
+- antiSMASH region GBKs commonly rebase extracted regions to coordinate zero. Such files do not automatically share their original chromosome coordinate system and must not use coordinate matching solely because they came from the same run.
```

Mirror the four modes and coordinate warning in [`SKILL.md`](SKILL.md), not just in the
reference file.

### Maintenance-command diff

The skill should use the checkout-safe pytest command already required by the plan:

```diff
diff --git a/SKILL.md b/SKILL.md
@@
-python -m pytest -q
-python -m pytest --cov=antismash_review --cov-report=term-missing
+python -m pytest -p no:cacheprovider --basetemp=.pytest_temp -q
+python -m pytest -p no:cacheprovider --basetemp=.pytest_temp \
+  --cov=antismash_review --cov-report=term-missing -q
```

## 6. Close the private-data acceptance tests

The manual audit confirmed the required behavior, but the automated private suite does not yet
assert the complete CLI attachment path. Add these assertions while keeping the directory
optional and untracked:

```diff
diff --git a/tests/test_sm_zpg19_integration.py b/tests/test_sm_zpg19_integration.py
@@
-from antismash_review.cli import main
+from antismash_review.cli import _load_review_records, main
@@
+def test_sm_zpg19_directory_loader_attaches_exactly_30_results(
+    sm_zpg19_dir: Path,
+) -> None:
+    manifest = discover(sm_zpg19_dir)
+    records, _ = _load_review_records(manifest, lenient=False)
+
+    assert len(records) == 1
+    assert len(records[0].clusterblast_results) == 30
+    assert all(result.source_format == "text" for result in records[0].clusterblast_results)
+    assert len(
+        {
+            (result.record_id, result.region_number, result.search_type)
+            for result in records[0].clusterblast_results
+        }
+    ) == 30
```

Prefer exposing a non-underscored shared loader if tests need to import it. Do not duplicate CLI
loading logic inside the test.

Also automate the already confirmed region-only case: parse the ten region GBKs, merge text and
JSON results, attach them, and assert three unique search types on each region record.

## 7. Release validation

The source checkout passes `pip check`, but the isolated wheel gate was not executable because
the host interpreter reports:

```text
ModuleNotFoundError: No module named 'setuptools'
WARNING: Package(s) not found: build, setuptools, wheel
```

This is an environment blocker, not proof of a packaging failure. Repeat the gate in a clean
environment containing the declared build backend:

```powershell
python -m venv .release-venv
.\.release-venv\Scripts\python -m pip install --upgrade pip setuptools wheel build
.\.release-venv\Scripts\python -m build
$wheel = Get-ChildItem dist\antismash_review-*.whl | Select-Object -First 1
.\.release-venv\Scripts\python -m pip install $wheel.FullName
.\.release-venv\Scripts\python -c "import antismash_review; print(antismash_review.__version__)"
.\.release-venv\Scripts\antismash-review --help
.\.release-venv\Scripts\python -m pip check
```

Use an external temporary directory or add `.release-venv/` to the local ignore configuration;
do not stage generated environments, wheels, or the private SM-ZPG19 bundle.

The package version `0.1.0` and record schema version `0.2.0` are separate version domains. Do
not change the package version merely to make it numerically equal to the schema version; make a
release-version decision explicitly after the defects above are fixed.

## Final quality gate

Run after all fixes:

```powershell
python -m ruff check .
python -m ruff format --check antismash_review tests
python -m mypy antismash_review
python -m pytest -p no:cacheprovider --basetemp=.pytest_temp -q
python -m pytest -p no:cacheprovider --basetemp=.pytest_temp `
  --cov=antismash_review --cov-report=term-missing -q
git -c safe.directory='D:/W/Skills Claude/antiSMASH-review' diff --check
```

Repeat all CLI smoke tests from the second-refactor plan. Final acceptance requires:

- no input representation can be overwritten;
- every supported sidecar failure is either fatal with status 2 or explicit as a lenient
  diagnostic;
- coordinate matching uses record identity independently of display IDs;
- coverage is at least 93%;
- documentation and output consistently say **record matching mode**;
- private integration tests pass when the data is present and skip cleanly otherwise;
- Markdown, build, import, entry-point, and `pip check` gates pass;
- the working tree contains no generated outputs or private biological inputs.
