from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from antismash_review.cli import main
from antismash_review.compare import (
    compare_records,
    intergenic_summary,
)
from antismash_review.exporters.compare_json import dumps_comparison
from antismash_review.exporters.compare_markdown import render_comparison
from antismash_review.genbank import parse_genbank
from antismash_review.models import (
    CollectionFeature,
    Gene,
    Location,
    LocationPart,
    Record,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _minimal_record(record_id: str = "REC1", **kwargs: object) -> Record:
    defaults: dict[str, object] = {
        "record_id": record_id,
        "name": record_id,
        "description": "test record",
        "length": 1000,
        "molecule_type": "DNA",
        "topology": "linear",
        "source_path": Path("/data/test.gbk"),
        "source_sha256": "abc123sha",
        "antismash_version": "8.0.4",
        "organism": "synthetic",
        "taxonomy": [],
    }
    defaults.update(kwargs)
    return Record(**defaults)  # type: ignore[arg-type]


def _simple_location(start: int, end: int) -> Location:
    part = LocationPart(start=start, end=end, strand=1)
    return Location(
        start=start,
        end=end,
        strand=1,
        parts=(part,),
        cross_origin=False,
        original=f"{start}..{end}",
    )


def _make_gene(locus_tag: str, start: int, end: int, kind: str = "biosynthetic") -> Gene:
    return Gene(
        location=_simple_location(start, end),
        locus_tag=locus_tag,
        gene=None,
        product=None,
        protein_id=None,
        translation=None,
        gene_kind=kind,
        gene_functions=[],
        ec_numbers=[],
        db_xrefs=[],
        notes=[],
        inference=[],
        region_numbers=[],
        candidate_cluster_numbers=[],
        protocluster_numbers=[],
        proto_core_numbers=[],
        qualifiers={},
    )


def _make_region(number: int, start: int, end: int, products: list[str]) -> CollectionFeature:
    return CollectionFeature(
        feature_type="region",
        number=number,
        location=_simple_location(start, end),
        products=products,
        references=[],
        kind=None,
        category=None,
        rules=[],
        smiles=[],
        polymer=[],
        core_location=None,
        cutoff=None,
        neighbourhood=None,
        creating_tool=None,
        contig_edge=None,
        qualifiers={},
    )


def test_compare_identical_records_zero_deltas() -> None:
    rec1 = _minimal_record("REC1")
    rec1.regions.append(_make_region(1, 0, 500, ["prodA", "prodB"]))
    rec1.genes.append(_make_gene("G1", 10, 100))

    rec2 = _minimal_record("REC1")
    rec2.regions.append(_make_region(1, 0, 500, ["prodA", "prodB"]))
    rec2.genes.append(_make_gene("G1", 10, 100))

    result = compare_records([rec1], [rec2], left_input=Path("l.gb"), right_input=Path("r.gb"))
    assert len(result.matched) == 1
    comp = result.matched[0]
    assert comp.left_gene_count == comp.right_gene_count == 1
    assert comp.left_region_count == comp.right_region_count == 1
    assert comp.gained_products == []
    assert comp.lost_products == []
    assert comp.new_diagnostics == []
    assert comp.resolved_diagnostics == []
    assert result.unmatched_left == []
    assert result.unmatched_right == []


def test_compare_product_multiplicity_and_order() -> None:
    rec1 = _minimal_record("REC1")
    rec1.regions.append(_make_region(1, 0, 500, ["prodA", "prodA", "prodB"]))

    rec2 = _minimal_record("REC1")
    rec2.regions.append(_make_region(1, 0, 500, ["prodA", "prodC", "prodC"]))

    result = compare_records([rec1], [rec2], left_input=Path("l.gb"), right_input=Path("r.gb"))
    comp = result.matched[0]
    assert comp.gained_products == ["prodC", "prodC"]
    assert comp.lost_products == ["prodA", "prodB"]


def test_compare_record_id_mode_duplicate_rejection() -> None:
    rec1a = _minimal_record("REC1")
    rec1b = _minimal_record("REC1")
    rec2 = _minimal_record("REC1")

    with pytest.raises(ValueError, match="Duplicate record ID found on left"):
        compare_records([rec1a, rec1b], [rec2], left_input=Path("l.gb"), right_input=Path("r.gb"))


def test_compare_record_region_mode() -> None:
    rec1 = _minimal_record("REC1")
    rec1.regions.append(_make_region(1, 0, 500, ["prodA"]))

    rec2 = _minimal_record("REC1")
    rec2.regions.append(_make_region(1, 0, 500, ["prodA"]))

    result = compare_records(
        [rec1],
        [rec2],
        left_input=Path("l.gb"),
        right_input=Path("r.gb"),
        match_method="record_region",
    )
    assert len(result.matched) == 1
    assert result.matched[0].match_key == "REC1:region_1"

    # Multi-region record rejected
    rec_multi = _minimal_record("REC1")
    rec_multi.regions.append(_make_region(1, 0, 200, ["p1"]))
    rec_multi.regions.append(_make_region(2, 300, 500, ["p2"]))
    with pytest.raises(ValueError, match="does not contain exactly one numbered region"):
        compare_records(
            [rec_multi],
            [rec2],
            left_input=Path("l.gb"),
            right_input=Path("r.gb"),
            match_method="record_region",
        )


def test_compare_single_record_mode() -> None:
    rec1 = _minimal_record("LEFT_ID")
    rec2 = _minimal_record("RIGHT_ID")

    result = compare_records(
        [rec1],
        [rec2],
        left_input=Path("l.gb"),
        right_input=Path("r.gb"),
        match_method="single_record",
    )
    assert len(result.matched) == 1
    assert result.matched[0].left_record_id == "LEFT_ID"
    assert result.matched[0].right_record_id == "RIGHT_ID"
    assert result.matched[0].match_key == "LEFT_ID <-> RIGHT_ID"

    # Multiple records on one side fails
    with pytest.raises(ValueError, match="single_record matching requires exactly one record"):
        compare_records(
            [rec1, _minimal_record("EXTRA")],
            [rec2],
            left_input=Path("l.gb"),
            right_input=Path("r.gb"),
            match_method="single_record",
        )


def test_compare_coordinate_overlap_requires_assumption_flag() -> None:
    rec1 = _minimal_record("REC1")
    rec2 = _minimal_record("REC2")

    with pytest.raises(ValueError, match="requires explicit --assume-shared-coordinate-system"):
        compare_records(
            [rec1],
            [rec2],
            left_input=Path("l.gb"),
            right_input=Path("r.gb"),
            match_method="coordinate_overlap",
            assume_shared_coordinate_system=False,
        )


def test_compare_coordinate_overlap_matching_and_evidence() -> None:
    rec1 = _minimal_record("LEFT")
    rec1.genes.append(_make_gene("G1", 0, 100))

    rec2 = _minimal_record("RIGHT")
    rec2.genes.append(_make_gene("G2", 0, 100))

    result = compare_records(
        [rec1],
        [rec2],
        left_input=Path("l.gb"),
        right_input=Path("r.gb"),
        match_method="coordinate_overlap",
        assume_shared_coordinate_system=True,
        min_reciprocal_overlap=0.80,
    )
    assert len(result.matched) == 1
    comp = result.matched[0]
    assert comp.coordinate_evidence is not None
    ev = comp.coordinate_evidence
    assert ev.overlap_bp == 100
    assert ev.left_span_bp == 100
    assert ev.right_span_bp == 100
    assert ev.left_overlap_fraction == 1.0
    assert ev.right_overlap_fraction == 1.0


def test_compare_coordinate_overlap_below_threshold_and_ambiguity() -> None:
    rec1 = _minimal_record("LEFT")
    rec1.genes.append(_make_gene("G1", 0, 100))

    rec2 = _minimal_record("RIGHT")
    rec2.genes.append(_make_gene("G2", 50, 150))  # overlap is 50/100 = 0.50

    result = compare_records(
        [rec1],
        [rec2],
        left_input=Path("l.gb"),
        right_input=Path("r.gb"),
        match_method="coordinate_overlap",
        assume_shared_coordinate_system=True,
        min_reciprocal_overlap=0.80,
    )
    assert len(result.matched) == 0
    assert result.unmatched_left == ["LEFT"]
    assert result.unmatched_right == ["RIGHT"]


def test_coordinate_overlap_allows_repeated_display_ids() -> None:
    left_a = _minimal_record("REPEATED")
    left_a.genes.append(_make_gene("LA", 0, 100))
    left_b = _minimal_record("REPEATED")
    left_b.genes.append(_make_gene("LB", 200, 300))

    right_a = _minimal_record("RIGHT")
    right_a.genes.append(_make_gene("RA", 0, 100))
    right_b = _minimal_record("RIGHT")
    right_b.genes.append(_make_gene("RB", 200, 300))

    result = compare_records(
        [left_a, left_b],
        [right_a, right_b],
        left_input=Path("left.gb"),
        right_input=Path("right.gb"),
        match_method="coordinate_overlap",
        assume_shared_coordinate_system=True,
    )

    assert len(result.matched) == 2
    assert result.unmatched_left == []
    assert result.unmatched_right == []
    evidence = [item.coordinate_evidence for item in result.matched]
    assert all(item is not None for item in evidence)
    assert [item.overlap_bp for item in evidence if item is not None] == [100, 100]


def test_coordinate_overlap_rejects_equal_candidate_tie() -> None:
    left = _minimal_record("LEFT")
    left.genes.append(_make_gene("L", 0, 100))
    right_a = _minimal_record("RIGHT_A")
    right_a.genes.append(_make_gene("RA", 0, 100))
    right_b = _minimal_record("RIGHT_B")
    right_b.genes.append(_make_gene("RB", 0, 100))

    with pytest.raises(ValueError, match="Ambiguous coordinate match"):
        compare_records(
            [left],
            [right_a, right_b],
            left_input=Path("left.gb"),
            right_input=Path("right.gb"),
            match_method="coordinate_overlap",
            assume_shared_coordinate_system=True,
        )


def test_coordinate_overlap_rejects_multiple_left_to_one_right() -> None:
    left_a = _minimal_record("LEFT_A")
    left_a.genes.append(_make_gene("LA", 0, 100))
    left_b = _minimal_record("LEFT_B")
    left_b.genes.append(_make_gene("LB", 0, 100))
    right = _minimal_record("RIGHT")
    right.genes.append(_make_gene("R", 0, 100))

    with pytest.raises(ValueError, match="Non-one-to-one coordinate match"):
        compare_records(
            [left_a, left_b],
            [right],
            left_input=Path("left.gb"),
            right_input=Path("right.gb"),
            match_method="coordinate_overlap",
            assume_shared_coordinate_system=True,
        )


@pytest.mark.parametrize("match_method", ["record_id", "record_region", "single_record"])
def test_assumption_flag_rejected_for_every_non_coordinate_mode(
    match_method: str,
) -> None:
    sem = FIXTURES / "semantics.gb"
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "compare",
                str(sem),
                str(sem),
                "--match-by",
                match_method,
                "--assume-shared-coordinate-system",
            ]
        )
    assert exc.value.code == 2


def test_intergenic_summary_handles_cross_origin_cds() -> None:
    record = parse_genbank(FIXTURES / "cross-origin.gb")[0]
    summary = intergenic_summary(record)
    assert summary.circular_wrap_included
    assert summary.gap_count == 1
    assert summary.total_bp == 60


def test_compare_record_id_right_duplicate_rejected() -> None:
    rec1 = _minimal_record("REC1")
    rec2a = _minimal_record("REC1")
    rec2b = _minimal_record("REC1")
    with pytest.raises(ValueError, match="Duplicate record ID found on right"):
        compare_records([rec1], [rec2a, rec2b], left_input=Path("l.gb"), right_input=Path("r.gb"))


def test_compare_single_record_empty_rejected() -> None:
    with pytest.raises(ValueError, match="single_record matching requires exactly one record"):
        compare_records(
            [],
            [_minimal_record("R")],
            left_input=Path("l.gb"),
            right_input=Path("r.gb"),
            match_method="single_record",
        )


def test_compare_record_region_missing_number_rejected() -> None:
    rec1 = _minimal_record("REC1")
    rec2 = _minimal_record("REC2")
    rec2.regions.append(_make_region(1, 0, 500, ["p"]))
    with pytest.raises(ValueError, match="does not contain exactly one numbered region"):
        compare_records(
            [rec1],
            [rec2],
            left_input=Path("l.gb"),
            right_input=Path("r.gb"),
            match_method="record_region",
        )


def test_compare_unknown_match_method_rejected() -> None:
    rec1 = _minimal_record("REC1")
    rec2 = _minimal_record("REC2")
    with pytest.raises(ValueError, match="Unknown match method"):
        compare_records(
            [rec1],
            [rec2],
            left_input=Path("l.gb"),
            right_input=Path("r.gb"),
            match_method="invalid_mode",
        )


def test_intergenic_summary_linear_and_circular() -> None:
    # Linear with two non-overlapping genes: (0, 10) and (20, 30) -> gap is 10 bp
    rec_lin = _minimal_record("LIN", length=100, topology="linear")
    rec_lin.genes.append(_make_gene("G1", 0, 10))
    rec_lin.genes.append(_make_gene("G2", 20, 30))
    summary_lin = intergenic_summary(rec_lin)
    assert summary_lin.gap_count == 1
    assert summary_lin.total_bp == 10
    assert summary_lin.mean_bp == 10.0
    assert summary_lin.median_bp == 10.0
    assert summary_lin.max_bp == 10
    assert not summary_lin.circular_wrap_included

    # Circular with same genes: wrap gap is (100 - 30) + 0 = 70 bp
    rec_circ = _minimal_record("CIRC", length=100, topology="circular")
    rec_circ.genes.append(_make_gene("G1", 0, 10))
    rec_circ.genes.append(_make_gene("G2", 20, 30))
    summary_circ = intergenic_summary(rec_circ)
    assert summary_circ.gap_count == 2
    assert summary_circ.total_bp == 80  # 10 + 70
    assert summary_circ.max_bp == 70
    assert summary_circ.circular_wrap_included


def test_comparison_json_and_markdown_exports() -> None:
    rec1 = parse_genbank(FIXTURES / "semantics.gb")[0]
    rec2 = parse_genbank(FIXTURES / "semantics.gb")[0]
    result = compare_records(
        [rec1], [rec2], left_input=FIXTURES / "semantics.gb", right_input=FIXTURES / "semantics.gb"
    )

    json_str = dumps_comparison(result)
    doc = json.loads(json_str)
    assert doc["schema_name"] == "antismash-review-comparison"
    assert doc["schema_version"] == "0.1.0"
    assert len(doc["comparison"]["matched"]) == 1

    md_str = render_comparison(result)
    assert "# antiSMASH comparative review" in md_str
    assert "SEMANTICS.1" in md_str
    assert "Regions | 1 | 1 | 0" in md_str


def test_cli_compare_command(capsys: object) -> None:
    sem = FIXTURES / "semantics.gb"
    assert main(["compare", str(sem), str(sem)]) == 0
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "# antiSMASH comparative review" in captured.out


def test_cli_compare_coordinate_validation(capsys: object) -> None:
    sem = FIXTURES / "semantics.gb"
    # coordinate_overlap without flag fails
    with pytest.raises(SystemExit) as exc1:
        main(["compare", str(sem), str(sem), "--match-by", "coordinate_overlap"])
    assert exc1.value.code == 2

    # flag without coordinate_overlap fails
    with pytest.raises(SystemExit) as exc2:
        main(
            [
                "compare",
                str(sem),
                str(sem),
                "--match-by",
                "record_id",
                "--assume-shared-coordinate-system",
            ]
        )
    assert exc2.value.code == 2

    # bad threshold fails (0.0)
    with pytest.raises(SystemExit) as exc3:
        main(
            [
                "compare",
                str(sem),
                str(sem),
                "--match-by",
                "coordinate_overlap",
                "--assume-shared-coordinate-system",
                "--min-reciprocal-overlap",
                "0.0",
            ]
        )
    assert exc3.value.code == 2

    # bad threshold fails (1.5)
    with pytest.raises(SystemExit) as exc4:
        main(
            [
                "compare",
                str(sem),
                str(sem),
                "--match-by",
                "coordinate_overlap",
                "--assume-shared-coordinate-system",
                "--min-reciprocal-overlap",
                "1.5",
            ]
        )
    assert exc4.value.code == 2


def test_cli_compare_overwrite_protection(tmp_path: Path, capsys: object) -> None:
    src_l = tmp_path / "left.gb"
    src_r = tmp_path / "right.gb"
    shutil.copy(FIXTURES / "semantics.gb", src_l)
    shutil.copy(FIXTURES / "semantics.gb", src_r)

    # Trying to overwrite either left or right input fails
    assert main(["compare", str(src_l), str(src_r), "--output", str(src_l)]) == 2
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "refusing to overwrite input file" in captured.err

    assert main(["compare", str(src_l), str(src_r), "--output", str(src_r)]) == 2
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "refusing to overwrite input file" in captured.err

    # Unselected region GenBank in a directory input is also protected
    left_dir = tmp_path / "left_dir"
    left_dir.mkdir()
    shutil.copy(FIXTURES / "semantics.gb", left_dir / "sample.gbk")
    region_file = left_dir / "sample.region001.gbk"
    shutil.copy(FIXTURES / "semantics.gb", region_file)
    original_region = region_file.read_bytes()

    assert main(["compare", str(left_dir), str(src_r), "--output", str(region_file)]) == 2
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "refusing to overwrite input file" in captured.err
    assert region_file.read_bytes() == original_region


def test_cli_compare_to_output_file(tmp_path: Path) -> None:
    sem = FIXTURES / "semantics.gb"
    out_file = tmp_path / "comparison.md"
    assert main(["compare", str(sem), str(sem), "--output", str(out_file)]) == 0
    assert out_file.exists()
    content = out_file.read_text(encoding="utf-8")
    assert "# antiSMASH comparative review" in content


def test_compare_markdown_full_features_and_unmatched() -> None:
    rec1 = _minimal_record("REC1")
    rec1.regions.append(_make_region(1, 0, 500, ["prodA", "prodB"]))
    rec1.genes.append(_make_gene("G1", 10, 100))

    rec2 = _minimal_record("REC1")
    rec2.regions.append(_make_region(1, 0, 500, ["prodA", "prodC"]))
    rec2.genes.append(_make_gene("G1", 10, 100))

    rec_unmatched_l = _minimal_record("EXTRA_LEFT")
    rec_unmatched_r = _minimal_record("EXTRA_RIGHT")

    result = compare_records(
        [rec1, rec_unmatched_l],
        [rec2, rec_unmatched_r],
        left_input=Path("left.gb"),
        right_input=Path("right.gb"),
        match_method="record_id",
    )
    assert result.unmatched_left == ["EXTRA_LEFT"]
    assert result.unmatched_right == ["EXTRA_RIGHT"]

    md = render_comparison(result)
    assert "- Gained products: prodC" in md
    assert "- Lost products: prodB" in md
    assert "## Unmatched left records" in md
    assert "- `EXTRA_LEFT`" in md
    assert "## Unmatched right records" in md
    assert "- `EXTRA_RIGHT`" in md


def test_compare_markdown_single_record_and_diagnostics() -> None:
    rec1 = _minimal_record("REC1")
    rec2 = _minimal_record("REC2")

    result = compare_records(
        [rec1],
        [rec2],
        left_input=Path("left.gb"),
        right_input=Path("right.gb"),
        match_method="single_record",
    )
    md = render_comparison(result)
    assert "- Note: single_record is an explicit user-requested pairing" in md


def test_compare_cli_format_json_and_errors(tmp_path: Path, capsys: object) -> None:
    sem = FIXTURES / "semantics.gb"

    # compare --format json
    assert main(["compare", str(sem), str(sem), "--format", "json"]) == 0
    captured = capsys.readouterr()  # type: ignore[union-attr]
    doc = json.loads(captured.out)
    assert doc["schema_name"] == "antismash-review-comparison"

    # Nonexistent file returns status 2
    assert main(["compare", "/nonexistent/left.gb", str(sem)]) == 2
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "error:" in captured.err

    # Corrupt file returns status 2
    bad = tmp_path / "bad.gb"
    bad.write_text("NOT A GENBANK", encoding="utf-8")
    assert main(["compare", str(bad), str(sem)]) == 2
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "error:" in captured.err


def test_compare_json_encoder_default_error() -> None:
    from antismash_review.exporters.compare_json import _json_default

    with pytest.raises(TypeError):
        _json_default(object())
