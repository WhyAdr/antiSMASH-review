# antiSMASH-review fourth-pass fixes

This plan contains the remaining work found after commit `c7c2710`. The three defects from
[`antiSMASH-review-3rd-fix.md`](antiSMASH-review-3rd-fix.md) are fixed: overwrite protection,
public structural ClusterBlast errors, and index-based coordinate matching all passed their
regressions. The current baseline is 112 passing tests with 94% coverage, eight passing private
SM-ZPG19 tests, and ten passing CLI smoke tests when pytest receives a fresh accessible
temporary base.

Do not broaden this work into full native antiSMASH JSON-to-`Record` parsing or cross-isolate
homology matching.

## Remaining work

| Priority | Area | Required outcome |
|---|---|---|
| 1 | ClusterBlast JSON scalar validation | Reject typed-schema violations instead of coercing, retaining, or discarding them silently |
| 2 | Repeatable pytest command | Stop reusing the inaccessible checkout-local `.pytest_temp` path |
| 3 | Documentation | Correct the optional overlap-threshold wording and mark the third-pass plan as implemented |
| 4 | Release gate | Build, install, import, invoke, and validate the wheel in a build-capable environment |

## 1. Close strict ClusterBlast scalar validation

### Observed residual behavior

The public structural error boundary works, but structurally valid JSON can still violate the
typed models without raising `ClusterBlastParseError`:

```text
missing_container_record_id -> accepted
region_number: true          -> retained as True
total_hits: "ten"            -> retained as a string
hits: "one"                  -> retained as a string
subject_index: "zero"        -> silently changed to None
```

This conflicts with the model types and the strict sidecar contract:

- `ClusterBlastResult.record_id: str`
- `ClusterBlastResult.region_number: int`
- `ClusterBlastResult.total_hits: int | None`
- integer score/count fields on `ClusterBlastHit`
- `ClusterBlastPairing.subject_index: int | None`

Python booleans require explicit handling because `isinstance(True, int)` is true.

### Add scalar validators

Add centralized validators near `ClusterBlastParseError` in
[`antismash_review/clusterblast.py`](antismash_review/clusterblast.py):

```diff
diff --git a/antismash_review/clusterblast.py b/antismash_review/clusterblast.py
@@
 import hashlib
 import json
+import math
 import re
@@
 class ClusterBlastParseError(RuntimeError):
     """A ClusterBlast sidecar was recognized but could not be parsed safely."""


+def _required_string(value: object, field: str, path: Path) -> str:
+    if not isinstance(value, str):
+        raise ClusterBlastParseError(
+            f"ClusterBlast {field} is missing or not a string in {path}"
+        )
+    return value
+
+
+def _required_nonempty_string(value: object, field: str, path: Path) -> str:
+    result = _required_string(value, field, path)
+    if not result:
+        raise ClusterBlastParseError(f"ClusterBlast {field} is empty in {path}")
+    return result
+
+
+def _optional_string(value: object, field: str, path: Path) -> str | None:
+    if value is None:
+        return None
+    return _required_string(value, field, path)
+
+
+def _required_integer(value: object, field: str, path: Path) -> int:
+    if type(value) is not int:
+        raise ClusterBlastParseError(
+            f"ClusterBlast {field} is not an integer in {path}: {value!r}"
+        )
+    return value
+
+
+def _optional_integer(value: object, field: str, path: Path) -> int | None:
+    if value is None:
+        return None
+    return _required_integer(value, field, path)
+
+
+def _required_float(value: object, field: str, path: Path) -> float:
+    if type(value) not in {int, float}:
+        raise ClusterBlastParseError(
+            f"ClusterBlast {field} is not numeric in {path}: {value!r}"
+        )
+    result = float(value)
+    if not math.isfinite(result):
+        raise ClusterBlastParseError(
+            f"ClusterBlast {field} is not finite in {path}: {value!r}"
+        )
+    return result
+
+
+def _optional_float(value: object, field: str, path: Path) -> float | None:
+    if value is None:
+        return None
+    return _required_float(value, field, path)
+
+
 def _sha256(path: Path) -> str:
```

Do not use permissive `int(value)`, `float(value)`, or `str(value)` conversions for JSON fields.
Those conversions hide schema drift. Text parsing may continue converting lexical numeric fields
because text input is necessarily string-valued.

### Validate record and result identity

```diff
diff --git a/antismash_review/clusterblast.py b/antismash_review/clusterblast.py
@@
-        rec_id = rec.get("id")
+        rec_id = _required_nonempty_string(
+            rec.get("id"), "containing record id", path
+        )
@@
-        cb_record_id = cb_mod.get("record_id")
-        if not isinstance(cb_record_id, str):
-            raise ClusterBlastParseError(
-                f"ClusterBlast module record_id is missing or not a string in {path}"
-            )
-        if rec_id is not None and cb_record_id != rec_id:
+        cb_record_id = _required_nonempty_string(
+            cb_mod.get("record_id"), "module record_id", path
+        )
+        if cb_record_id != rec_id:
             raise ClusterBlastParseError(
                 f"ClusterBlast module record_id {cb_record_id!r} does not match "
                 f"record id {rec_id!r} in {path}"
             )
@@
-                region_number = region_res["region_number"]
-                if not isinstance(region_number, int):
+                region_number = _required_integer(
+                    region_res["region_number"], "region_number", path
+                )
+                if region_number < 1:
                     raise ClusterBlastParseError(
-                        f"ClusterBlast region_number is not an integer in {path}"
+                        f"ClusterBlast region_number must be positive in {path}: "
+                        f"{region_number}"
                     )
-                total_hits = region_res.get("total_hits")
+                total_hits = _optional_integer(
+                    region_res.get("total_hits"), "total_hits", path
+                )
+                if total_hits is not None and total_hits < 0:
+                    raise ClusterBlastParseError(
+                        f"ClusterBlast total_hits must not be negative in {path}: "
+                        f"{total_hits}"
+                    )
```

Requiring a containing record ID implements the existing requirement to validate module
`record_id` against the containing JSON record ID. Absence is not a successful comparison.

### Validate pairing fields

```diff
diff --git a/antismash_review/clusterblast.py b/antismash_review/clusterblast.py
@@
                         parts = query_str.split("|", 5)
-                        if len(parts) < 5:
+                        if len(parts) < 5 or not parts[4]:
                             raise ClusterBlastParseError(
                                 f"Malformed query string in {path}: {query_str!r}"
                             )
                         query_gene = parts[4]
@@
                         pairings.append(
                             ClusterBlastPairing(
                                 query_gene=query_gene,
-                                subject_gene=str(pairing_dict["name"]),
-                                percent_identity=float(pairing_dict["perc_ident"]),
-                                blast_score=float(pairing_dict["blastscore"]),
-                                percent_coverage=float(pairing_dict["perc_coverage"]),
-                                evalue=float(pairing_dict["evalue"]),
-                                subject_protein_id=pairing_dict.get("locus_tag"),
-                                subject_index=subject_idx if isinstance(subject_idx, int) else None,
+                                subject_gene=_required_nonempty_string(
+                                    pairing_dict.get("name"), "pairing name", path
+                                ),
+                                percent_identity=_required_float(
+                                    pairing_dict.get("perc_ident"), "pairing perc_ident", path
+                                ),
+                                blast_score=_required_float(
+                                    pairing_dict.get("blastscore"), "pairing blastscore", path
+                                ),
+                                percent_coverage=_required_float(
+                                    pairing_dict.get("perc_coverage"),
+                                    "pairing perc_coverage",
+                                    path,
+                                ),
+                                evalue=_required_float(
+                                    pairing_dict.get("evalue"), "pairing evalue", path
+                                ),
+                                subject_protein_id=_optional_string(
+                                    pairing_dict.get("locus_tag"),
+                                    "pairing locus_tag",
+                                    path,
+                                ),
+                                subject_index=_optional_integer(
+                                    subject_idx, "pairing subject_index", path
+                                ),
                             )
                         )
```

If the supported antiSMASH schema requires `subject_index` rather than allowing null, use
`_required_integer` instead. Confirm this against the local antiSMASH 8.0.4 JSON before choosing;
do not infer it from the Python model alone.

### Validate hit fields

```diff
diff --git a/antismash_review/clusterblast.py b/antismash_review/clusterblast.py
@@
-                    raw_blast_score = hit_details.get("blast_score")
                     rankings.append(
                         ClusterBlastHit(
                             rank=rank_idx,
-                            accession=str(hit_info["accession"]),
-                            description=str(hit_info["description"]),
-                            cluster_type=hit_info.get("cluster_type"),
-                            num_hits=hit_details.get("hits"),
-                            core_gene_hits=hit_details.get("core_gene_hits"),
-                            blast_score=float(raw_blast_score)
-                            if raw_blast_score is not None
-                            else None,
-                            synteny_score=hit_details.get("synteny_score"),
-                            core_bonus=hit_details.get("core_bonus"),
-                            similarity=hit_details.get("similarity"),
+                            accession=_required_nonempty_string(
+                                hit_info.get("accession"), "hit accession", path
+                            ),
+                            description=_required_string(
+                                hit_info.get("description"), "hit description", path
+                            ),
+                            cluster_type=_optional_string(
+                                hit_info.get("cluster_type"), "hit cluster_type", path
+                            ),
+                            num_hits=_optional_integer(
+                                hit_details.get("hits"), "hit count", path
+                            ),
+                            core_gene_hits=_optional_integer(
+                                hit_details.get("core_gene_hits"),
+                                "core gene hit count",
+                                path,
+                            ),
+                            blast_score=_optional_float(
+                                hit_details.get("blast_score"), "hit blast_score", path
+                            ),
+                            synteny_score=_optional_integer(
+                                hit_details.get("synteny_score"), "synteny_score", path
+                            ),
+                            core_bonus=_optional_integer(
+                                hit_details.get("core_bonus"), "core_bonus", path
+                            ),
+                            similarity=_optional_integer(
+                                hit_details.get("similarity"), "similarity", path
+                            ),
                             pairings=pairings,
                         )
                     )
```

Before enforcing non-negative ranges for every score, inspect the supported schema and real
fixture. Counts and similarity should not be negative, but score semantics should not be guessed.

### Regression tests

Add a small valid-document factory to [`tests/test_clusterblast.py`](tests/test_clusterblast.py)
and mutate one scalar per case:

```diff
diff --git a/tests/test_clusterblast.py b/tests/test_clusterblast.py
@@
+def _minimal_clusterblast_json_document() -> dict[str, object]:
+    return {
+        "records": [
+            {
+                "id": "contig_1",
+                "modules": {
+                    "antismash.modules.clusterblast": {
+                        "schema_version": 2,
+                        "record_id": "contig_1",
+                        "general": {
+                            "schema_version": 5,
+                            "results": [
+                                {
+                                    "region_number": 1,
+                                    "total_hits": 1,
+                                    "ranking": [
+                                        [
+                                            {
+                                                "accession": "ACC1",
+                                                "description": "description",
+                                                "cluster_type": "NRPS",
+                                            },
+                                            {
+                                                "hits": 1,
+                                                "core_gene_hits": 0,
+                                                "blast_score": 5.0,
+                                                "synteny_score": 1,
+                                                "core_bonus": 0,
+                                                "similarity": 50,
+                                                "pairings": [
+                                                    [
+                                                        "input|c1|1-100|+|QUERY|NRPS",
+                                                        0,
+                                                        {
+                                                            "name": "subject",
+                                                            "perc_ident": 90.0,
+                                                            "blastscore": 100.0,
+                                                            "perc_coverage": 95.0,
+                                                            "evalue": 1e-20,
+                                                            "locus_tag": "SUBJECT_1",
+                                                        },
+                                                    ]
+                                                ],
+                                            },
+                                        ]
+                                    ],
+                                }
+                            ],
+                        },
+                    }
+                },
+            }
+        ]
+    }
```

Use focused tests rather than one opaque mega-parameterization:

```diff
diff --git a/tests/test_clusterblast.py b/tests/test_clusterblast.py
@@
+def test_json_clusterblast_requires_containing_record_id(tmp_path: Path) -> None:
+    document = _minimal_clusterblast_json_document()
+    del document["records"][0]["id"]  # type: ignore[index]
+    path = tmp_path / "missing_record_id.json"
+    path.write_text(json.dumps(document), encoding="utf-8")
+
+    with pytest.raises(ClusterBlastParseError, match="containing record id"):
+        parse_clusterblast_json(path)
+
+
+@pytest.mark.parametrize("value", [True, False, "1", 1.0, None])
+def test_json_clusterblast_rejects_non_integer_region_number(
+    tmp_path: Path,
+    value: object,
+) -> None:
+    document = _minimal_clusterblast_json_document()
+    result = document["records"][0]["modules"][  # type: ignore[index]
+        "antismash.modules.clusterblast"
+    ]["general"]["results"][0]
+    result["region_number"] = value
+    path = tmp_path / "bad_region.json"
+    path.write_text(json.dumps(document), encoding="utf-8")
+
+    with pytest.raises(ClusterBlastParseError, match="region_number"):
+        parse_clusterblast_json(path)
+
+
+@pytest.mark.parametrize(
+    ("field", "value"),
+    [
+        ("total_hits", "ten"),
+        ("hits", "one"),
+        ("core_gene_hits", True),
+        ("synteny_score", 1.5),
+        ("core_bonus", "zero"),
+        ("similarity", False),
+    ],
+)
+def test_json_clusterblast_rejects_non_integer_counters(
+    tmp_path: Path,
+    field: str,
+    value: object,
+) -> None:
+    document = _minimal_clusterblast_json_document()
+    result = document["records"][0]["modules"][  # type: ignore[index]
+        "antismash.modules.clusterblast"
+    ]["general"]["results"][0]
+    if field == "total_hits":
+        result[field] = value
+    else:
+        result["ranking"][0][1][field] = value
+    path = tmp_path / f"bad_{field}.json"
+    path.write_text(json.dumps(document), encoding="utf-8")
+
+    with pytest.raises(ClusterBlastParseError):
+        parse_clusterblast_json(path)
+
+
+def test_json_clusterblast_rejects_subject_index_coercion(tmp_path: Path) -> None:
+    document = _minimal_clusterblast_json_document()
+    pairing = document["records"][0]["modules"][  # type: ignore[index]
+        "antismash.modules.clusterblast"
+    ]["general"]["results"][0]["ranking"][0][1]["pairings"][0]
+    pairing[1] = "zero"
+    path = tmp_path / "bad_subject_index.json"
+    path.write_text(json.dumps(document), encoding="utf-8")
+
+    with pytest.raises(ClusterBlastParseError, match="subject_index"):
+        parse_clusterblast_json(path)
```

Add parallel tests for required strings, optional strings, non-finite numeric values, and valid
`None` values. Re-run the real SM-ZPG19 JSON to ensure the stricter checks match antiSMASH 8.0.4
rather than assumptions.

### Acceptance

- Every accepted JSON value conforms to its dataclass annotation.
- Booleans are rejected for integer and floating-point fields.
- Invalid non-null values are never changed to `None`.
- Missing containing record IDs and mismatched module IDs fail publicly.
- Strict CLI mode returns status 2 without a traceback.
- Lenient CLI mode emits `clusterblast_parse_failed`.
- All 30 SM-ZPG19 JSON results still parse with the established totals and parity.

## 2. Make pytest basetemp repeatable

### Finding

The fixed `.pytest_temp` path is currently inaccessible in this checkout. The documented command
produced 38 setup errors with `PermissionError: [WinError 5]`, while fresh basetemp directories
produced 112 passing tests and 94% coverage. Even `Get-Acl .pytest_temp` was denied.

This is path/ownership state, not a package-test failure. The maintenance command should not
depend on reusing a directory that may have been created by a differently privileged process.

### SKILL.md diff

Use a per-process system temporary directory on Windows:

````diff
diff --git a/SKILL.md b/SKILL.md
@@
-```bash
+```powershell
 python -m ruff check .
 python -m ruff format --check antismash_review tests
 python -m mypy antismash_review
-python -m pytest -p no:cacheprovider --basetemp=.pytest_temp -q
-python -m pytest -p no:cacheprovider --basetemp=.pytest_temp --cov=antismash_review --cov-report=term-missing -q
+$pytestBase = Join-Path $env:TEMP "antismash-review-pytest-$PID"
+python -m pytest -p no:cacheprovider --basetemp=$pytestBase -q
+python -m pytest -p no:cacheprovider --basetemp=$pytestBase `
+  --cov=antismash_review --cov-report=term-missing -q
 python -m antismash_review --help
 ```
````

For cross-platform documentation, give separate PowerShell and POSIX examples rather than
embedding shell-specific environment syntax in one block. A random UUID is even safer than a
PID if multiple runs may share the same shell process.

Do not attempt to take ownership of or delete the current inaccessible `.pytest_temp` directory
from package code. Repository code must not manage checkout ACLs.

### Acceptance

- The documented commands pass twice consecutively from the same checkout.
- They do not create untracked files inside the repository.
- A stale or inaccessible old basetemp does not affect a new run.
- Test failures are distinguishable from temporary-directory setup failures.

## 3. Correct documentation status and option wording

### Coordinate threshold wording

`--assume-shared-coordinate-system` is mandatory in coordinate mode.
`--min-reciprocal-overlap` is optional and defaults to 0.80.

```diff
diff --git a/SKILL.md b/SKILL.md
@@
-   Supported record matching modes: `record_id` (default), `record_region`, `single_record`, and `coordinate_overlap` (requires `--assume-shared-coordinate-system` and `--min-reciprocal-overlap`). Use coordinate matching only when coordinate correspondence is independently established (e.g. re-annotations of identical contigs), not between arbitrary isolates or rebased region files.
+   Supported record matching modes: `record_id` (default), `record_region`, `single_record`, and `coordinate_overlap`. Coordinate mode requires `--assume-shared-coordinate-system`; `--min-reciprocal-overlap` optionally changes the default 0.80 threshold. Use coordinate matching only when coordinate correspondence is independently established (e.g. re-annotations of identical contigs), not between arbitrary isolates or rebased region files.
```

### Mark the third-pass plan as implemented

```diff
diff --git a/antiSMASH-review-3rd-fix.md b/antiSMASH-review-3rd-fix.md
@@
-The snippets below are implementation-ready guidance. They have not been applied to the
-package code.
+The snippets below were implementation guidance for commit `c7c2710`. The three principal
+behavioral fixes were applied and verified; remaining scalar-validation and release work is
+tracked in [`4th-fix.md`](4th-fix.md).
```

Keep the third-pass document as an audit record; do not rewrite its observed evidence or make it
look as though the original defects never existed.

## 4. Complete the isolated release gate

The current interpreter passes `pip check` but lacks `setuptools`, `wheel`, and `build`, so the
wheel gate remains unverified. Run this in a disposable environment where dependency installation
is authorized:

```powershell
python -m venv .release-venv
$releasePython = (Resolve-Path .release-venv\Scripts\python.exe).Path
& $releasePython -m pip install --upgrade pip setuptools wheel build
& $releasePython -m build
$wheel = Get-ChildItem dist\antismash_review-*.whl | Select-Object -First 1
& $releasePython -m pip install $wheel.FullName
$entryPoint = (Resolve-Path .release-venv\Scripts\antismash-review.exe).Path
$outsideCheckout = Join-Path $env:TEMP "antismash-review-wheel-check-$PID"
New-Item -ItemType Directory -Force $outsideCheckout | Out-Null
Push-Location $outsideCheckout
try {
  & $releasePython -c "import antismash_review; print(antismash_review.__version__)"
  & $entryPoint --help
  & $releasePython -m pip check
} finally {
  Pop-Location
}
```

Validate both the import and installed console entry point from outside the source checkout so
the working tree cannot mask missing wheel contents. Confirm that `SKILL.md` and references are
distributed by the intended skill-delivery mechanism; they are not Python package modules.

Generated `.release-venv/`, `dist/`, `build/`, and `*.egg-info/` paths must remain ignored and
untracked.

## Final verification

Run the updated repeatable quality gate, then:

```powershell
python -m antismash_review --help
python -m antismash_review inspect tests/fixtures/semantics.gb
python -m antismash_review inspect tests/fixtures/semantics.gb --format json
python -m antismash_review inspect tests/fixtures/semantics.gb --format tsv
python -m antismash_review inspect tests/fixtures/semantics.gb --format gene-tsv
python -m antismash_review inspect tests/fixtures/semantics.gb --format domain-tsv
python -m antismash_review inspect SM-ZPG19--NOACCESSION-antismash --format clusterblast-tsv
python -m antismash_review compare tests/fixtures/semantics.gb `
  tests/fixtures/semantics.gb --match-by record_id
python -m antismash_review compare tests/fixtures/semantics.gb `
  tests/fixtures/semantics.gb --match-by coordinate_overlap `
  --assume-shared-coordinate-system
python -m antismash_review compare --help
```

Final acceptance requires:

- full suite passes with coverage at least 93%;
- all private SM-ZPG19 tests pass when the bundle is available and skip cleanly otherwise;
- every accepted ClusterBlast JSON scalar conforms to the public dataclass type;
- every malformed supported sidecar becomes `ClusterBlastParseError` or an explicit lenient
  diagnostic;
- both documented pytest invocations work repeatedly without repository-local ACL dependence;
- coordinate documentation distinguishes the required assumption flag from the optional
  threshold override;
- Ruff, formatting, mypy, Markdown checks, `git diff --check`, wheel build/install/import,
  console entry point, and `pip check` all pass;
- the working tree contains no generated outputs or private biological inputs.
