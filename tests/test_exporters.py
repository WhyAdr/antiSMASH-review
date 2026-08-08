"""Tests for exporter edge cases: contig_edge tri-state, Markdown/TSV escaping."""

from pathlib import Path

from antismash_review.exporters.markdown import render_records
from antismash_review.exporters.tables import render_tsv
from antismash_review.genbank import parse_genbank
from antismash_review.models import (
    CollectionFeature,
    Location,
    LocationPart,
    Record,
)

ROOT = Path(__file__).resolve().parent / "fixtures"


def _minimal_record(**kwargs: object) -> Record:
    defaults: dict[str, object] = {
        "record_id": "TEST.1",
        "name": "TEST",
        "description": "test record",
        "length": 100,
        "molecule_type": "DNA",
        "topology": "linear",
        "source_path": Path("/synthetic/test.gb"),
        "source_sha256": "0" * 64,
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


def test_contig_edge_true() -> None:
    record = _minimal_record()
    record.regions.append(
        CollectionFeature(
            feature_type="region",
            number=1,
            location=_simple_location(0, 100),
            products=["test"],
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
            contig_edge=True,
            qualifiers={},
        )
    )
    tsv = render_tsv([record])
    lines = tsv.strip().splitlines()
    assert lines[1].split("\t")[3] == "true"


def test_contig_edge_false() -> None:
    record = _minimal_record()
    record.regions.append(
        CollectionFeature(
            feature_type="region",
            number=1,
            location=_simple_location(0, 100),
            products=["test"],
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
    tsv = render_tsv([record])
    lines = tsv.strip().splitlines()
    assert lines[1].split("\t")[3] == "false"


def test_contig_edge_none_is_empty() -> None:
    record = _minimal_record()
    record.regions.append(
        CollectionFeature(
            feature_type="region",
            number=1,
            location=_simple_location(0, 100),
            products=["test"],
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
    )
    tsv = render_tsv([record])
    lines = tsv.strip().splitlines()
    assert lines[1].split("\t")[3] == ""


def test_contig_edge_mixed_true_wins() -> None:
    record = _minimal_record()
    for edge in [False, True]:
        record.regions.append(
            CollectionFeature(
                feature_type="region",
                number=1,
                location=_simple_location(0, 100),
                products=["test"],
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
                contig_edge=edge,
                qualifiers={},
            )
        )
    tsv = render_tsv([record])
    lines = tsv.strip().splitlines()
    assert lines[1].split("\t")[3] == "true"


def test_markdown_escaping() -> None:
    """Products with special characters are rendered safely in Markdown."""
    record = parse_genbank(ROOT / "semantics.gb")[0]
    md = render_records([record])
    assert "# antiSMASH review" in md
    assert "hybrid-a" in md
    assert "hybrid-b" in md


def test_tsv_header_order() -> None:
    """TSV header columns are in the documented order."""
    record = parse_genbank(ROOT / "semantics.gb")[0]
    tsv = render_tsv([record])
    header = tsv.splitlines()[0]
    expected_fields = [
        "filename",
        "record_id",
        "region_products",
        "contig_edge",
        "core_genes",
        "total_genes",
        "nrps_pks_domains",
        "all_domains",
        "modules",
        "raw_pfam_hits",
        "deduplicated_pfam_hits",
        "diagnostics",
    ]
    assert header.split("\t") == expected_fields


def test_tsv_row_cell_count() -> None:
    """Every TSV data row has the same number of cells as the header."""
    record = parse_genbank(ROOT / "semantics.gb")[0]
    tsv = render_tsv([record])
    lines = tsv.strip().splitlines()
    header_count = len(lines[0].split("\t"))
    for line in lines[1:]:
        assert len(line.split("\t")) == header_count
