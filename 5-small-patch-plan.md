# antiSMASH-review fifth small-patch plan

This plan addresses the GLM audit of the current `antiSMASH-review` checkout. It is a
plan-only artifact: the diffs below are implementation guidance and are not applied by
this patch.

## Audit baseline

The audit was checked against commit `af23c85` (`feat: expose typed Python API`) and the
working tree's current source.

| Check | Current evidence | Consequence |
|---|---|---|
| Tests | `16 passed` | The passing count is mostly API/schema smoke coverage. |
| Total coverage | `39%` | `cli.py`, `discovery.py`, and `compare.py` are at `0%`; `genbank.py` is at `41%`; `clusterblast.py` is at `41%`. |
| Public fixture | No committed GenBank fixture exists | Core parser behavior is not portable across clones. |
| Private integration data | `SM-ZPG19--NOACCESSION-antismash/` is ignored | Keep it untracked; do not make it the only regression authority. |
| Sidecar layout | The local tree has 10 top-level `.txt` files in each of `clusterblast/`, `knownclusterblast/`, and `subclusterblast/`; nested `knownclusterblast/regionN/` directories contain HTML | Preserve non-recursive `.txt` discovery and make that contract executable. |
| Coordinate mode | Default threshold is `0.80`; coordinate matching already uses right-record positions internally | Add threshold and repeated-ID tests; do not change the default based on an unverified isolate comparison. |
| Specificity | Local GBKs contain repeated `KR activity: inactive` and `KR stereochemistry: C2` values | The current `list[str]` plus raw qualifier map retains these values; add a regression test instead of introducing an unverified T2PKS schema. |
| Dead alias | `rg` finds only `_load_review_records = load_review_records` in `cli.py` | Remove the alias after confirming no import depends on it. |

## Scope and acceptance target

The patch should make the portable checkout useful on a fresh clone without importing any
unpublished biological input. It should:

- add a source-verified synthetic Biopython fixture builder and tests for GenBank parsing,
  review diagnostics, discovery, comparison, ClusterBlast enrichment, exports, and CLI
  dispatch;
- emit a high-value `pseudogene_in_cluster` warning for a standalone `/gene` feature with
  `/pseudo` that overlaps a `region`;
- emit an aggregated `unrecognized_feature_type` notice for non-structural feature types
  retained only in `raw_features`;
- codify the current non-recursive sidecar layout rather than guessing at antiSMASH output
  variants that are not present in the verified tree;
- prove the `0.80` reciprocal-overlap boundary, rejection below the boundary, ambiguity
  behavior, and the explicit shared-coordinate-system guard;
- remove the unused underscored loader alias; and
- prove that repeated domain specificity qualifiers remain present in both the model and
  exported text.

The coverage gate should be raised from the current `39%` to at least `80%` overall, with
`genbank.py`, `discovery.py`, `compare.py`, and `cli.py` each at least `80%`. No module named
in the audit may remain at `0%`. A later release patch can tighten the total percentage after the
full exporter and private integration matrix is restored.

## Patch 1 - Add a portable synthetic fixture and core parser tests

Use only synthetic sequence and qualifier values. The fixture must be built by Biopython and
must exercise the qualifier names already verified in real antiSMASH GenBank files:

- one `region`, `cand_cluster`, `protocluster`, and `proto_core` hierarchy;
- one fuzzy edge `CDS` and one standalone pseudo `gene`;
- one NRPS/PKS `aSDomain`, one duplicate domain ID, and repeated specificity values;
- one `aSModule` with one valid and one dangling domain reference;
- one `CDS_motif` without an e-value;
- two version-drifted `PFAM_domain` hits that deduplicate;
- one `tRNA`, one `rRNA`, and one `repeat_region` retained as unknown raw features; and
- a required `source` feature, which is structural and must not create an unknown-feature
  notice by itself.

### Fixture builder diff

```diff
diff --git a/tests/__init__.py b/tests/__init__.py
new file mode 100644
--- /dev/null
+++ b/tests/__init__.py
@@
+"""Tests for the antiSMASH review package."""

diff --git a/tests/fixtures/__init__.py b/tests/fixtures/__init__.py
new file mode 100644
--- /dev/null
+++ b/tests/fixtures/__init__.py
@@
+"""Synthetic, non-biological test fixtures."""

diff --git a/tests/fixtures/build_fixture.py b/tests/fixtures/build_fixture.py
new file mode 100644
--- /dev/null
+++ b/tests/fixtures/build_fixture.py
@@
+from __future__ import annotations
+
+from pathlib import Path
+
+from Bio import SeqIO
+from Bio.Seq import Seq
+from Bio.SeqFeature import BeforePosition, ExactPosition, FeatureLocation, SeqFeature
+from Bio.SeqRecord import SeqRecord
+
+
+def _feature(
+    feature_type: str,
+    start: int,
+    end: int,
+    qualifiers: dict[str, list[str]] | None = None,
+    *,
+    fuzzy_start: bool = False,
+) -> SeqFeature:
+    left = BeforePosition(start) if fuzzy_start else ExactPosition(start)
+    location = FeatureLocation(left, ExactPosition(end), strand=1)
+    return SeqFeature(location=location, type=feature_type, qualifiers=qualifiers or {})
+
+
+def write_synthetic_genbank(path: Path) -> Path:
+    """Write a source-verified, biologically synthetic antiSMASH-like record."""
+    record = SeqRecord(
+        Seq("A" * 400),
+        id="SYNTH.1",
+        name="SYNTH.1",
+        description="synthetic antiSMASH review fixture",
+    )
+    record.annotations["molecule_type"] = "DNA"
+    record.annotations["topology"] = "linear"
+    record.annotations["structured_comment"] = {
+        "antiSMASH-Data": {"Version": "8.0.4"}
+    }
+    record.features = [
+        _feature("source", 0, 400, {"organism": ["synthetic bacterium"]}),
+        _feature("region", 0, 400, {
+            "region_number": ["1"],
+            "candidate_cluster_numbers": ["1", "1"],
+            "product": ["NRPS", "synthetic product"],
+            "contig_edge": ["true"],
+        }),
+        _feature("cand_cluster", 20, 380, {
+            "candidate_cluster_number": ["1"],
+            "protoclusters": ["1"],
+            "product": ["NRPS"],
+        }),
+        _feature("protocluster", 20, 380, {
+            "protocluster_number": ["1"],
+            "candidate_cluster_numbers": ["1"],
+        }),
+        _feature("proto_core", 0, 120, {"protocluster_number": ["1"]}),
+        _feature("CDS", 0, 90, {
+            "locus_tag": ["SYN_CDS_1"],
+            "gene_kind": ["biosynthetic"],
+            "translation": ["M" * 30],
+        }, fuzzy_start=True),
+        _feature("gene", 100, 160, {
+            "locus_tag": ["SYN_PSEUDO_1"],
+            "pseudo": [""],
+            "note": ["frameshift introduced in synthetic fixture"],
+        }),
+        _feature("aSDomain", 30, 80, {
+            "domain_id": ["D1"],
+            "aSDomain": ["A"],
+            "aSTool": ["nrps_pks_domains"],
+            "locus_tag": ["SYN_CDS_1"],
+            "specificity": ["KR activity: inactive", "KR stereochemistry: C2"],
+        }),
+        _feature("aSDomain", 80, 120, {
+            "domain_id": ["D1"],
+            "aSDomain": ["TIGRFAM-domain"],
+            "aSTool": ["tigrfam"],
+            "locus_tag": ["SYN_CDS_1"],
+        }),
+        _feature("aSModule", 30, 120, {
+            "domains": ["D1", "MISSING"],
+            "locus_tags": ["SYN_CDS_1", "MISSING_CDS"],
+            "complete": [""],
+        }),
+        _feature("CDS_motif", 130, 150, {
+            "locus_tag": ["SYN_CDS_1"],
+            "core_sequence": ["SYNCORE"],
+        }),
+        _feature("PFAM_domain", 30, 80, {
+            "db_xref": ["PF00001.1"],
+            "locus_tag": ["SYN_CDS_1"],
+            "protein_start": ["1"],
+            "protein_end": ["10"],
+        }),
+        _feature("PFAM_domain", 30, 80, {
+            "db_xref": ["PF00001.2"],
+            "locus_tag": ["SYN_CDS_1"],
+            "protein_start": ["1"],
+            "protein_end": ["10"],
+        }),
+        _feature("tRNA", 170, 190),
+        _feature("rRNA", 200, 220),
+        _feature("repeat_region", 230, 250),
+    ]
+    SeqIO.write(record, path, "genbank")
+    return path
```

The builder must remain test-only. Do not copy a private `.gbk` into `tests/fixtures`, do not
embed a real isolate sequence, and do not add the ignored SM-ZPG19 directory to Git.

### Parser/review/export regression diff

```diff
diff --git a/tests/test_genbank_review.py b/tests/test_genbank_review.py
new file mode 100644
--- /dev/null
+++ b/tests/test_genbank_review.py
@@
+from pathlib import Path
+
+from antismash_review.exporters.entity_tables import render_domain_tsv
+from antismash_review.exporters.json_export import dumps_records
+from antismash_review.exporters.markdown import render_records
+from antismash_review.exporters.tables import render_tsv
+from antismash_review.genbank import parse_genbank
+from antismash_review.review import review_record
+from tests.fixtures.build_fixture import write_synthetic_genbank
+
+
+def test_synthetic_fixture_preserves_hierarchy_raw_features_and_qualifiers(
+    tmp_path: Path,
+) -> None:
+    record = parse_genbank(write_synthetic_genbank(tmp_path / "synthetic.gbk"))[0]
+
+    assert len(record.regions) == 1
+    assert len(record.candidate_clusters) == 1
+    assert len(record.protoclusters) == 1
+    assert len(record.proto_cores) == 1
+    assert len(record.genes) == 1  # standalone gene is raw evidence, not a CDS
+    assert [domain.tool for domain in record.domains] == [
+        "nrps_pks_domains",
+        "tigrfam",
+    ]
+    assert record.domains[0].specificity == [
+        "KR activity: inactive",
+        "KR stereochemistry: C2",
+    ]
+    assert record.domains[0].qualifiers["specificity"] == tuple(record.domains[0].specificity)
+    assert len(record.modules[0].missing_domain_ids) == 1
+    assert len(record.deduplicated_pfam_hits) == 1
+    assert {raw.feature_type for raw in record.raw_features} >= {
+        "gene",
+        "tRNA",
+        "rRNA",
+        "repeat_region",
+    }
+
+    codes = {diagnostic.code for diagnostic in review_record(record)}
+    assert "pseudogene_in_cluster" in codes
+    assert "unrecognized_feature_type" in codes
+    assert "partial_cds_at_edge" in codes
+    assert "module_domain_missing" in codes
+    assert "orphan_module_locus" in codes
+
+    assert "pseudogene_in_cluster" in dumps_records([record])
+    assert "pseudogene_in_cluster" in render_records([record])
+    assert "pseudogene_in_cluster" in render_tsv([record])
+    assert "KR stereochemistry: C2" in render_domain_tsv([record])
```

The two package marker files above make the fixture import deterministic across pytest
invocations and installed environments.

## Patch 2 - Add evidence-scoped parser and review diagnostics

Do not promote every `gene` feature into the CDS-backed `Gene` model. A standalone gene can
lack a translation and is semantically different from a CDS. Keep it in `raw_features`, then
review it by coordinate overlap.

### Adapter set and unknown-feature diagnostic diff

```diff
diff --git a/antismash_review/genbank.py b/antismash_review/genbank.py
@@
 _GENE_FUNCTION_RE = re.compile(
     r"^(?P<category>[\w-]+)(?:\s+\((?P<tool>[^)]+)\))?(?:\s+(?P<description>.*))?$"
 )
+
+_ADAPTED_FEATURE_TYPES = frozenset(
+    {
+        "region",
+        "cand_cluster",
+        "protocluster",
+        "proto_core",
+        "CDS",
+        "aSDomain",
+        "aSModule",
+        "CDS_motif",
+        "PFAM_domain",
+    }
+)
+_STRUCTURAL_RAW_ONLY_FEATURE_TYPES = frozenset({"source"})
@@
 def _adapt_feature(record: Record, raw: RawFeature) -> None:
@@
     elif raw.feature_type == "PFAM_domain":
         record.pfam_hits.append(_pfam_hit(raw))
+
+
+def _diagnose_unrecognized_features(record: Record) -> None:
+    unknown = sorted(
+        {
+            raw.feature_type
+            for raw in record.raw_features
+            if raw.feature_type not in _ADAPTED_FEATURE_TYPES
+            and raw.feature_type not in _STRUCTURAL_RAW_ONLY_FEATURE_TYPES
+        }
+    )
+    if unknown:
+        record.diagnostics.append(
+            Diagnostic(
+                code="unrecognized_feature_type",
+                severity=Severity.NOTICE,
+                message=(
+                    "Feature types retained only as raw evidence and not adapted: "
+                    + ", ".join(unknown)
+                ),
+                source=str(record.source_path),
+                record_id=record.record_id,
+            )
+        )
@@
     _resolve_modules(record)
     _assign_gene_memberships(record)
+    _diagnose_unrecognized_features(record)
     return record
```

The set must be kept next to the adapter dispatch so a new adapter cannot silently omit the
completeness check. The diagnostic is aggregated by feature type, not emitted once per raw
feature, so a genome with hundreds of repeated unknown features remains readable.

### Pseudo-gene review diff

```diff
diff --git a/antismash_review/review.py b/antismash_review/review.py
@@
-from .models import Diagnostic, Record, Severity
+from .models import Diagnostic, Location, RawFeature, Record, Severity
@@
+def _locations_overlap(left: Location, right: Location) -> bool:
+    return any(
+        left_part.start < right_part.end and right_part.start < left_part.end
+        for left_part in left.parts
+        for right_part in right.parts
+    )
+
+
+def _raw_label(raw: RawFeature) -> str:
+    for key in ("locus_tag", "gene"):
+        values = raw.qualifiers.get(key, ())
+        if values and values[0]:
+            return values[0]
+    return "unlabelled gene"
+
+
 def _extend_consistency_diagnostics(
     record: Record,
     diagnostics: list[Diagnostic],
 ) -> None:
+    for raw in record.raw_features:
+        if raw.feature_type != "gene" or "pseudo" not in raw.qualifiers:
+            continue
+        overlapping_regions = [
+            region
+            for region in record.regions
+            if _locations_overlap(raw.location, region.location)
+        ]
+        if not overlapping_regions:
+            continue
+        region_numbers = sorted(
+            region.number for region in overlapping_regions if region.number is not None
+        )
+        region_text = ", ".join(str(number) for number in region_numbers) or "unnumbered"
+        diagnostics.append(
+            Diagnostic(
+                code="pseudogene_in_cluster",
+                severity=Severity.WARNING,
+                message=(
+                    f"Pseudo gene {_raw_label(raw)} overlaps antiSMASH region(s): "
+                    f"{region_text}; inspect frameshift or annotation evidence"
+                ),
+                source=str(record.source_path),
+                record_id=record.record_id,
+                feature_index=raw.feature_index,
+            )
+        )
+
     gene_locus_tags = {gene.locus_tag for gene in record.genes if gene.locus_tag}
```

The `/pseudo` qualifier is detected by key presence because GenBank writers may serialize a
flag as an empty value, `true`, or another textual value. The diagnostic says only that a
pseudo gene overlaps a region; it must not claim a proven frameshift, lost function, or
compound identity.

### Contract documentation diff

```diff
diff --git a/references/semantic-contract.md b/references/semantic-contract.md
@@
 ## Review boundaries
@@
 - Do not emit a generic motif-confidence warning based on one e-value cutoff.
+- Standalone `gene` features remain raw evidence. A `/pseudo` gene overlapping a region
+  emits `pseudogene_in_cluster`; this does not prove a frameshift or functional loss.
+- Feature types outside the adapter set are retained in `raw_features` and produce one
+  aggregated `unrecognized_feature_type` NOTICE per record, except the structural `source`
+  feature.
@@
 | `missing_nrps_pks_architecture` | WARNING | Region products contain NRPS/PKS terms but no aSTool=nrps_pks_domains domains were parsed | Does not prove the annotation is wrong; the relevant domains may be in a different record |
+| `pseudogene_in_cluster` | WARNING | A standalone `/gene` with `/pseudo` overlaps one or more regions | Does not prove a frameshift, functional loss, or compound identity |
+| `unrecognized_feature_type` | NOTICE | A non-structural feature type is retained only in `raw_features` | Does not mean the source annotation is invalid; it marks adapter coverage |
```

Mirror the user-facing diagnostic names in `SKILL.md` under the conservative interpretation
section, and state that `gene-tsv` remains one row per CDS-backed `Gene`; raw-only features
are visible through JSON/Markdown diagnostics and the raw feature model.

## Patch 3 - Exercise discovery, sidecars, CLI, and comparison

### Resolve the sidecar recursion question with a fixture test

The verified local result tree settles the question for this checkout: recognized text files
are direct children of the three canonical sidecar directories. Nested `knownclusterblast`
folders contain HTML detail assets, not parser inputs. Do not switch `glob("*.txt")` to
`rglob("*.txt")` in this small patch. A future antiSMASH layout change should be a separate
compatibility decision with a real fixture.

```diff
diff --git a/tests/test_discovery.py b/tests/test_discovery.py
new file mode 100644
--- /dev/null
+++ b/tests/test_discovery.py
@@
+from pathlib import Path
+
+from antismash_review.discovery import discover
+
+
+def test_discover_classifies_genbanks_and_prefers_aggregate(tmp_path: Path) -> None:
+    (tmp_path / "aggregate.gbk").write_text("", encoding="utf-8")
+    (tmp_path / "contig.region001.gbk").write_text("", encoding="utf-8")
+    manifest = discover(tmp_path)
+
+    assert manifest.aggregate_genbanks == (tmp_path / "aggregate.gbk",)
+    assert manifest.region_genbanks == (tmp_path / "contig.region001.gbk",)
+
+
+def test_discover_sidecars_only_reads_canonical_top_level_text(tmp_path: Path) -> None:
+    for name in ("clusterblast", "knownclusterblast", "subclusterblast"):
+        sidecar = tmp_path / name
+        sidecar.mkdir()
+        direct = sidecar / "contig_1_c1.txt"
+        direct.write_text("direct", encoding="utf-8")
+        nested = sidecar / "region1"
+        nested.mkdir()
+        (nested / "details.txt").write_text("html-sidecar placeholder", encoding="utf-8")
+
+    manifest = discover(tmp_path)
+
+    assert manifest.clusterblast_files == (tmp_path / "clusterblast" / "contig_1_c1.txt",)
+    assert manifest.knownclusterblast_files == (
+        tmp_path / "knownclusterblast" / "contig_1_c1.txt",
+    )
+    assert manifest.subclusterblast_files == (
+        tmp_path / "subclusterblast" / "contig_1_c1.txt",
+    )
+    assert all("details.txt" not in str(path) for path in manifest.knownclusterblast_files)
```

Also add a short sentence to the discovery section of `references/semantic-contract.md`:

```diff
diff --git a/references/semantic-contract.md b/references/semantic-contract.md
@@
 - Parse text sidecars from canonical `clusterblast/`, `knownclusterblast/`, and `subclusterblast/` directories using natural region sorting (`_cN.txt`).
+- The current contract reads recognized `.txt` sidecars directly inside those canonical
+  directories; nested HTML/detail assets are ignored.
```

### Valid text sidecar and JSON enrichment tests

Add one tiny text fixture with synthetic identifiers and one valid JSON factory. Exercise
parsing, source provenance, precedence, empty-result retention, and attachment by
`(record_id, region_number, search_type)`. The text fixture may use this source-verified
shape:

```diff
diff --git a/tests/fixtures/clusterblast/contig_1_c1.txt b/tests/fixtures/clusterblast/contig_1_c1.txt
new file mode 100644
--- /dev/null
+++ b/tests/fixtures/clusterblast/contig_1_c1.txt
@@
+ClusterBlast scores for SYNTH.1
+Significant hits:
+1. SYNTH-HIT-1	Synthetic hit
+Details:
+1. SYNTH-HIT-1
+Type: NRPS
+Number of proteins with BLAST hits to this cluster: 1
+Cumulative BLAST score: 12.5
+Table of Blast hits
+SYN_CDS_1	SYNTH_SUBJECT	90.0	12.5	80.0	1e-10
+>>
```

The corresponding test diff should cover the public parser and merge/attach behavior:

```diff
diff --git a/tests/test_clusterblast_validation.py b/tests/test_clusterblast_validation.py
--- a/tests/test_clusterblast_validation.py
+++ b/tests/test_clusterblast_validation.py
@@
-from antismash_review.clusterblast import ClusterBlastParseError, parse_clusterblast_json
+from antismash_review.clusterblast import (
+    ClusterBlastParseError,
+    attach_clusterblast_results,
+    merge_clusterblast_results,
+    parse_clusterblast_text,
+)
+from antismash_review.genbank import parse_genbank
+from tests.fixtures.build_fixture import write_synthetic_genbank
+
+
+FIXTURE = Path(__file__).parent / "fixtures" / "clusterblast" / "contig_1_c1.txt"
+
+
+def test_text_clusterblast_parses_pairings_and_provenance() -> None:
+    result = parse_clusterblast_text(FIXTURE, search_type="clusterblast")
+
+    assert result.record_id == "SYNTH.1"
+    assert result.region_number == 1
+    assert result.source_format == "text"
+    assert result.rankings[0].pairings[0].query_gene == "SYN_CDS_1"
+    assert result.rankings[0].pairings[0].percent_coverage == 80.0
+
+
+def test_clusterblast_results_merge_text_precedence_and_attach(tmp_path: Path) -> None:
+    records = parse_genbank(write_synthetic_genbank(tmp_path / "synthetic.gbk"))
+    text_result = parse_clusterblast_text(FIXTURE, search_type="clusterblast")
+    merged = merge_clusterblast_results([text_result], [])
+    attach_clusterblast_results(records, merged)
+
+    assert len(records[0].clusterblast_results) == 1
+    assert records[0].clusterblast_results[0].source_format == "text"
+    assert records[0].clusterblast_results[0].source_sha256 == text_result.source_sha256
```

Keep the existing schema-validation tests. Extend them with valid zero-result JSON, all three
search sections, record/module mismatch, malformed container shapes, and one lenient CLI
diagnostic test. This is the coverage that proves the sidecar path rather than only its scalar
guard rails.

### CLI dispatch and output-format test diff

```diff
diff --git a/tests/test_cli.py b/tests/test_cli.py
new file mode 100644
--- /dev/null
+++ b/tests/test_cli.py
@@
+from pathlib import Path
+
+import pytest
+
+from antismash_review.cli import main
+from tests.fixtures.build_fixture import write_synthetic_genbank
+
+
+@pytest.mark.parametrize(
+    "output_format",
+    ["markdown", "json", "tsv", "gene-tsv", "domain-tsv", "clusterblast-tsv"],
+)
+def test_cli_inspect_dispatches_every_public_format(
+    tmp_path: Path,
+    capsys: pytest.CaptureFixture[str],
+    output_format: str,
+) -> None:
+    input_path = write_synthetic_genbank(tmp_path / "synthetic.gbk")
+
+    assert main(["inspect", str(input_path), "--format", output_format]) == 0
+    output = capsys.readouterr().out
+    assert output
+    if output_format == "json":
+        assert '"records"' in output
+    elif output_format == "markdown":
+        assert "# antiSMASH review" in output
+    else:
+        assert output.splitlines()[0]
+
+
+def test_cli_coordinate_mode_requires_explicit_shared_coordinates() -> None:
+    with pytest.raises(SystemExit) as exc_info:
+        main(["compare", "left.gbk", "right.gbk", "--match-by", "coordinate_overlap"])
+    assert exc_info.value.code == 2
```

Add focused CLI tests for malformed input status `2`, `--lenient` status `0` with an explicit
diagnostic, output-file writing, and refusing to overwrite an input. Use the already-public
`load_review_records`; do not import a private alias from tests.

### Coordinate threshold and repeated-display-ID test diff

The default remains `0.80`. Add tests that make the policy observable without using the
unavailable BK71-I annotation vintages. A real pair of re-annotations can be added later as
an optional private integration test, but it must not be required for the distributable suite.

```diff
diff --git a/tests/test_compare.py b/tests/test_compare.py
new file mode 100644
--- /dev/null
+++ b/tests/test_compare.py
@@
+from pathlib import Path
+
+import pytest
+
+from antismash_review.compare import compare_records
+from antismash_review.models import Gene, Location, LocationPart, Record
+
+
+def _record(record_id: str, start: int, end: int) -> Record:
+    location = Location(
+        start=start,
+        end=end,
+        strand=1,
+        parts=(LocationPart(start, end, 1),),
+        cross_origin=False,
+        original=f"{start}..{end}",
+    )
+    gene = Gene(
+        location=location,
+        locus_tag=f"{record_id}_{start}",
+        gene=None,
+        product=None,
+        protein_id=None,
+        translation=None,
+        gene_kind="unclassified",
+        gene_functions=[],
+        ec_numbers=[],
+        db_xrefs=[],
+        notes=[],
+        inference=[],
+        region_numbers=[],
+        candidate_cluster_numbers=[],
+        protocluster_numbers=[],
+        proto_core_numbers=[],
+        qualifiers={},
+    )
+    return Record(
+        record_id=record_id,
+        name=record_id,
+        description="synthetic comparison record",
+        length=1000,
+        molecule_type="DNA",
+        topology="linear",
+        source_path=Path(f"{record_id}.gbk"),
+        source_sha256="",
+        antismash_version=None,
+        organism=None,
+        taxonomy=[],
+        genes=[gene],
+    )
+
+
+def test_coordinate_overlap_matches_at_default_reciprocal_threshold() -> None:
+    result = compare_records(
+        [_record("LEFT", 0, 100)],
+        [_record("RIGHT", 20, 120)],
+        left_input=Path("left.gbk"),
+        right_input=Path("right.gbk"),
+        match_method="coordinate_overlap",
+        assume_shared_coordinate_system=True,
+    )
+
+    assert result.min_reciprocal_overlap == 0.80
+    assert result.matched[0].coordinate_evidence is not None
+    assert result.matched[0].coordinate_evidence.left_overlap_fraction == 0.80
+
+
+def test_coordinate_overlap_below_threshold_is_unmatched() -> None:
+    result = compare_records(
+        [_record("LEFT", 0, 100)],
+        [_record("RIGHT", 21, 121)],
+        left_input=Path("left.gbk"),
+        right_input=Path("right.gbk"),
+        match_method="coordinate_overlap",
+        assume_shared_coordinate_system=True,
+    )
+
+    assert result.matched == []
+    assert result.unmatched_left == ["LEFT"]
+    assert result.unmatched_right == ["RIGHT"]
+
+
+def test_coordinate_overlap_uses_record_positions_not_repeated_ids() -> None:
+    result = compare_records(
+        [_record("REPEATED", 0, 100), _record("REPEATED", 200, 300)],
+        [_record("REPEATED", 0, 100), _record("REPEATED", 200, 300)],
+        left_input=Path("left.gbk"),
+        right_input=Path("right.gbk"),
+        match_method="coordinate_overlap",
+        assume_shared_coordinate_system=True,
+    )
+
+    assert len(result.matched) == 2
+    assert result.unmatched_left == []
+    assert result.unmatched_right == []
+
+
+def test_coordinate_overlap_rejects_ambiguous_candidates() -> None:
+    with pytest.raises(ValueError, match="Ambiguous coordinate match"):
+        compare_records(
+            [_record("LEFT", 0, 100)],
+            [_record("RIGHT_A", 0, 100), _record("RIGHT_B", 0, 100)],
+            left_input=Path("left.gbk"),
+            right_input=Path("right.gbk"),
+            match_method="coordinate_overlap",
+            assume_shared_coordinate_system=True,
+        )
```

Add the remaining matching-mode matrix (`record_id`, `record_region`, `single_record`),
duplicate-ID failures where those modes require unique keys, circular intergenic summaries,
diagnostic/product deltas, and the true multiple-left-to-one-right failure. These tests should
assert both values and deterministic ordering of unmatched records.

### Remove the dead alias diff

```diff
diff --git a/antismash_review/cli.py b/antismash_review/cli.py
@@
-_load_review_records = load_review_records
```

Verify the cleanup with `rg -n "_load_review_records" .` and update any test that still imports
the old name to import `load_review_records`.

## Patch 4 - Lock specificity semantics and release documentation

The current `Domain` model does not discard the local T2PKS-adjacent evidence that was
verified in the real bundle: repeated `/specificity` qualifiers are parsed into an ordered
`list[str]`, and the original repeated values remain in `Domain.qualifiers`. The domain TSV
joins the list without dropping entries. No local fixture contains a distinct Minowa field,
so this patch must not invent one or claim that a Minowa prediction is parsed.

Add an explicit contract sentence and keep the future extension path clear:

```diff
diff --git a/references/semantic-contract.md b/references/semantic-contract.md
@@
 - Parse plural `domain_subtypes` and legacy singular `domain_subtype`.
+- Preserve repeated `/specificity` values in order in both `Domain.specificity` and the raw
+  qualifier map. The current model does not classify T2PKS/KR/Minowa semantics; downstream
+  code must interpret the retained strings conservatively.
```

```diff
diff --git a/SKILL.md b/SKILL.md
@@
 - Do not apply a universal motif e-value threshold. Interpret motif evidence using the producing tool and motif family.
+- Treat domain specificity as source evidence: repeated values, including KR activity or
+  stereochemistry strings, are retained but are not independently validated predictions.
```

If future data demonstrates that a structured Minowa field is being lost, add a typed parsed
view while retaining the existing raw list for backward compatibility. That is a separate
schema patch, not part of this small patch.

## Verification sequence

Run the following from a clean disposable pytest base directory:

```powershell
python -m ruff check .
python -m ruff format --check antismash_review tests
python -m mypy antismash_review
$planPytestBase = Join-Path $env:TEMP ("antismash-review-plan-" + [guid]::NewGuid())
python -m pytest -p no:cacheprovider --basetemp=$planPytestBase -q
python -m pytest -p no:cacheprovider --basetemp=$planPytestBase `
  --cov=antismash_review --cov-report=term-missing -q
python -m antismash_review --help
git -c safe.directory=* diff --check
```

Also run CLI smoke checks for all six inspect formats and the four comparison modes. If the
optional private bundle is present, run a separate ignored integration test that asserts the
known local facts (one aggregate record, 30 text sidecar results, and three search types per
region where applicable). If it is absent, the test must skip cleanly. Do not lower the public
coverage gate because private data is unavailable.

## Final acceptance

- `tests/fixtures/build_fixture.py` contains no real isolate sequence or private identifiers.
- Standalone pseudo genes overlapping regions produce `pseudogene_in_cluster` WARNINGs with
  source and raw-feature index evidence.
- Unrecognized non-structural feature types produce one `unrecognized_feature_type` NOTICE;
  `source` alone does not create noise.
- JSON, Markdown, and summary TSV expose diagnostics; domain TSV retains repeated specificity
  values.
- Sidecar discovery remains direct-child `.txt` discovery and is covered by a nested-file test.
- Coordinate mode remains explicit and defaults to `0.80`; exact-boundary, below-boundary,
  ambiguity, repeated-ID, and guard-flag tests pass.
- `_load_review_records` has no remaining references.
- Overall coverage is at least `80%`, and `genbank.py`, `discovery.py`, `compare.py`, and
  `cli.py` are each at least `80%`.
- Ruff, formatting, mypy, pytest, CLI help, and `git diff --check` pass.
- Private biological inputs and generated outputs remain untracked.
