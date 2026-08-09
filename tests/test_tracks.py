from __future__ import annotations

from pathlib import Path

from antismash_review.exporters.bed import render_bed
from antismash_review.exporters.gff3 import render_gff3
from antismash_review.genbank import parse_genbank
from antismash_review.locations import overlaps
from antismash_review.models import Location, LocationPart
from antismash_review.review import review_findings
from tests.fixtures.build_fixture import write_synthetic_genbank
from tests.test_assemblyline import _domain, _module, _record


def test_gff3_and_bed_keep_their_coordinate_contracts() -> None:
    record = _record(
        [_module(0, 100, ["CDS1"], ["D1"], strand=-1)],
        [_domain("D1", "AMP-binding", "CDS1")],
    )

    gff = render_gff3([record])
    bed = render_bed([record])

    assert "\taSModule\t1\t100\t.\t-\t.\t" in gff
    assert "ASSEMBLY.1\t0\t100\taSModule|" in bed


def test_rebased_region_seqid_and_duplicate_preferred_ids_are_deterministic() -> None:
    record = _record(
        [_module(0, 100, ["CDS1"], ["D1", "D1"])],
        [_domain("D1", "AMP-binding", "CDS1"), _domain("D1", "PCP", "CDS1")],
    )
    record.source_path = Path("contig_1.region003.gbk")

    gff = render_gff3([record])

    assert "contig_1.region003:ASSEMBLY.1\t" in gff
    assert "ID=ASSEMBLY.1:aSDomain:D1" in gff
    assert "ID=ASSEMBLY.1:aSDomain:D1.2" in gff


def test_compound_cross_origin_parts_and_attribute_escaping_are_preserved() -> None:
    record = _record(
        [_module(0, 100, ["CDS1"], ["D1"], pairings=["A;=, % -> X"])],
        [_domain("D1", "AMP-binding", "CDS1")],
    )
    compound = Location(
        start=0,
        end=100,
        strand=1,
        parts=(LocationPart(0, 20, 1), LocationPart(80, 100, 1)),
        cross_origin=True,
        original="join(1..20,81..100)",
    )
    record.modules[0].location = compound

    gff = render_gff3([record])

    assert "\taSModule\t1\t20\t" in gff
    assert "\taSModule\t81\t100\t" in gff
    assert "part=1%2F2" in gff
    assert "cross_origin=true" in gff
    assert "%3B%3D%2C%20%25" in gff


def test_localized_partial_cds_finding_is_exported_at_cds_location(tmp_path: Path) -> None:
    record = parse_genbank(write_synthetic_genbank(tmp_path / "synthetic.gbk"))[0]
    findings = review_findings(record)
    partial = next(
        finding for finding in findings if finding.diagnostic.code == "partial_cds_at_edge"
    )

    assert partial.location is not None
    assert partial.entity_type == "CDS"
    gff = render_gff3([record])
    assert "review_finding" in gff
    assert "code=partial_cds_at_edge" in gff
    assert "\treview_finding\t1\t90\t" in gff


def test_track_rendering_is_byte_deterministic_and_location_helper_is_half_open() -> None:
    record = _record(
        [_module(0, 100, ["CDS1"], ["D1"])],
        [_domain("D1", "AMP-binding", "CDS1")],
    )

    other = Location(
        100,
        200,
        1,
        (LocationPart(100, 200, 1),),
        False,
        "100..200",
    )
    assert overlaps(record.modules[0].location, other) is False
    assert render_gff3([record]) == render_gff3([record])
    assert render_bed([record]) == render_bed([record])
