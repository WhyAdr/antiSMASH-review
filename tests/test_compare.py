from __future__ import annotations

from pathlib import Path

import pytest

from antismash_review.compare import compare_records, intergenic_summary
from antismash_review.exporters.compare_json import dumps_comparison
from antismash_review.exporters.compare_markdown import render_comparison
from antismash_review.models import CollectionFeature, Gene, Location, LocationPart, Record


def _record(
    record_id: str,
    start: int,
    end: int,
    *,
    topology: str = "linear",
    region_number: int | None = None,
    products: list[str] | None = None,
) -> Record:
    location = Location(
        start=start,
        end=end,
        strand=1,
        parts=(LocationPart(start, end, 1),),
        cross_origin=False,
        original=f"{start}..{end}",
    )
    gene = Gene(
        location=location,
        locus_tag=f"{record_id}_{start}",
        gene=None,
        product=None,
        protein_id=None,
        translation=None,
        gene_kind="unclassified",
        gene_functions=[],
        ec_numbers=[],
        db_xrefs=[],
        notes=[],
        inference=[],
        region_numbers=[region_number] if region_number is not None else [],
        candidate_cluster_numbers=[],
        protocluster_numbers=[],
        proto_core_numbers=[],
        qualifiers={},
    )
    regions: list[CollectionFeature] = []
    if region_number is not None:
        regions.append(
            CollectionFeature(
                feature_type="region",
                number=region_number,
                location=location,
                products=products or ["NRPS"],
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
                contig_edge=False,
                qualifiers={},
            )
        )
    return Record(
        record_id=record_id,
        name=record_id,
        description="synthetic comparison record",
        length=1000,
        molecule_type="DNA",
        topology=topology,
        source_path=Path(f"{record_id}.gbk"),
        source_sha256="",
        antismash_version=None,
        organism=None,
        taxonomy=[],
        genes=[gene],
        regions=regions,
    )


def test_coordinate_overlap_matches_at_default_reciprocal_threshold() -> None:
    result = compare_records(
        [_record("LEFT", 0, 100)],
        [_record("RIGHT", 20, 120)],
        left_input=Path("left.gbk"),
        right_input=Path("right.gbk"),
        match_method="coordinate_overlap",
        assume_shared_coordinate_system=True,
    )

    assert result.min_reciprocal_overlap == 0.80
    assert len(result.matched) == 1
    assert result.matched[0].coordinate_evidence is not None
    assert result.matched[0].coordinate_evidence.left_overlap_fraction == 0.80
    assert result.matched[0].coordinate_evidence.right_overlap_fraction == 0.80


def test_coordinate_overlap_below_threshold_is_unmatched() -> None:
    result = compare_records(
        [_record("LEFT", 0, 100)],
        [_record("RIGHT", 21, 121)],
        left_input=Path("left.gbk"),
        right_input=Path("right.gbk"),
        match_method="coordinate_overlap",
        assume_shared_coordinate_system=True,
    )

    assert result.matched == []
    assert result.unmatched_left == ["LEFT"]
    assert result.unmatched_right == ["RIGHT"]


def test_coordinate_overlap_uses_record_positions_not_repeated_ids() -> None:
    result = compare_records(
        [_record("REPEATED", 0, 100), _record("REPEATED", 200, 300)],
        [_record("REPEATED", 0, 100), _record("REPEATED", 200, 300)],
        left_input=Path("left.gbk"),
        right_input=Path("right.gbk"),
        match_method="coordinate_overlap",
        assume_shared_coordinate_system=True,
    )

    assert len(result.matched) == 2
    assert result.unmatched_left == []
    assert result.unmatched_right == []


def test_coordinate_overlap_rejects_ambiguous_candidates() -> None:
    with pytest.raises(ValueError, match="Ambiguous coordinate match"):
        compare_records(
            [_record("LEFT", 0, 100)],
            [_record("RIGHT_A", 0, 100), _record("RIGHT_B", 0, 100)],
            left_input=Path("left.gbk"),
            right_input=Path("right.gbk"),
            match_method="coordinate_overlap",
            assume_shared_coordinate_system=True,
        )


def test_coordinate_overlap_rejects_non_one_to_one() -> None:
    with pytest.raises(ValueError, match="Non-one-to-one coordinate match"):
        compare_records(
            [_record("LEFT_A", 0, 100), _record("LEFT_B", 0, 100)],
            [_record("RIGHT", 0, 100)],
            left_input=Path("left.gbk"),
            right_input=Path("right.gbk"),
            match_method="coordinate_overlap",
            assume_shared_coordinate_system=True,
        )


def test_coordinate_overlap_requires_explicit_shared_flag() -> None:
    with pytest.raises(ValueError, match="requires explicit --assume-shared-coordinate-system"):
        compare_records(
            [_record("LEFT", 0, 100)],
            [_record("RIGHT", 0, 100)],
            left_input=Path("left.gbk"),
            right_input=Path("right.gbk"),
            match_method="coordinate_overlap",
            assume_shared_coordinate_system=False,
        )


def test_coordinate_overlap_invalid_threshold() -> None:
    with pytest.raises(ValueError, match="min_reciprocal_overlap must be in the interval"):
        compare_records(
            [_record("LEFT", 0, 100)],
            [_record("RIGHT", 0, 100)],
            left_input=Path("left.gbk"),
            right_input=Path("right.gbk"),
            match_method="coordinate_overlap",
            assume_shared_coordinate_system=True,
            min_reciprocal_overlap=0.0,
        )


def test_compare_record_id_matching() -> None:
    left = [_record("REC1", 0, 100), _record("REC2", 0, 100)]
    right = [_record("REC2", 0, 100), _record("REC3", 0, 100)]

    result = compare_records(
        left,
        right,
        left_input=Path("left.gbk"),
        right_input=Path("right.gbk"),
        match_method="record_id",
    )
    assert len(result.matched) == 1
    assert result.matched[0].match_key == "REC2"
    assert result.unmatched_left == ["REC1"]
    assert result.unmatched_right == ["REC3"]

    # Duplicate IDs on left
    with pytest.raises(ValueError, match="Duplicate record ID found on left"):
        compare_records(
            [_record("DUP", 0, 100), _record("DUP", 100, 200)],
            [_record("OTHER", 0, 100)],
            left_input=Path("left.gbk"),
            right_input=Path("right.gbk"),
            match_method="record_id",
        )

    # Duplicate IDs on right
    with pytest.raises(ValueError, match="Duplicate record ID found on right"):
        compare_records(
            [_record("OTHER", 0, 100)],
            [_record("DUP", 0, 100), _record("DUP", 100, 200)],
            left_input=Path("left.gbk"),
            right_input=Path("right.gbk"),
            match_method="record_id",
        )


def test_compare_record_region_matching() -> None:
    left = [_record("REC1", 0, 100, region_number=1), _record("REC1", 100, 200, region_number=2)]
    right = [_record("REC1", 0, 100, region_number=1), _record("REC1", 200, 300, region_number=3)]

    result = compare_records(
        left,
        right,
        left_input=Path("left.gbk"),
        right_input=Path("right.gbk"),
        match_method="record_region",
    )
    assert len(result.matched) == 1
    assert result.matched[0].match_key == "REC1:region_1"
    assert result.unmatched_left == ["REC1:region_2"]
    assert result.unmatched_right == ["REC1:region_3"]

    # Record without exactly one numbered region
    no_region = _record("NO_REG", 0, 100)
    with pytest.raises(ValueError, match="does not contain exactly one numbered region"):
        compare_records(
            [no_region],
            right,
            left_input=Path("left.gbk"),
            right_input=Path("right.gbk"),
            match_method="record_region",
        )


def test_compare_single_record_matching() -> None:
    left = [_record("LEFT_ID", 0, 100)]
    right = [_record("RIGHT_ID", 0, 100)]

    result = compare_records(
        left,
        right,
        left_input=Path("left.gbk"),
        right_input=Path("right.gbk"),
        match_method="single_record",
    )
    assert len(result.matched) == 1
    assert result.matched[0].match_key == "LEFT_ID <-> RIGHT_ID"

    # More than 1 record
    with pytest.raises(ValueError, match="single_record matching requires exactly one record"):
        compare_records(
            [_record("L1", 0, 100), _record("L2", 0, 100)],
            right,
            left_input=Path("left.gbk"),
            right_input=Path("right.gbk"),
            match_method="single_record",
        )


def test_compare_unknown_match_method() -> None:
    with pytest.raises(ValueError, match="Unknown match method"):
        compare_records(
            [_record("L", 0, 100)],
            [_record("R", 0, 100)],
            left_input=Path("left.gbk"),
            right_input=Path("right.gbk"),
            match_method="invalid_method",
        )


def test_intergenic_summary_circular_wrap() -> None:
    rec_circular = _record("CIRC", 100, 200, topology="circular")
    # length is 1000, gene is 100..200. Circular wrap gap is (1000 - 200) + 100 = 900
    summary = intergenic_summary(rec_circular)
    assert summary.circular_wrap_included is True
    assert summary.gap_count == 1
    assert summary.total_bp == 900


def test_compare_exporters() -> None:
    left = [_record("REC1", 0, 100, region_number=1, products=["NRPS", "T1PKS"])]
    right = [_record("REC1", 0, 100, region_number=1, products=["NRPS", "Terpene"])]

    result = compare_records(
        left,
        right,
        left_input=Path("left.gbk"),
        right_input=Path("right.gbk"),
        match_method="record_id",
    )
    assert result.matched[0].gained_products == ["Terpene"]
    assert result.matched[0].lost_products == ["T1PKS"]

    md = render_comparison(result)
    assert "# antiSMASH comparative review" in md
    assert "REC1" in md
    assert "Terpene" in md

    json_str = dumps_comparison(result)
    assert "antismash-review-comparison" in json_str
    assert "Terpene" in json_str
