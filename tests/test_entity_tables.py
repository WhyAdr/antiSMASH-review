from __future__ import annotations

import csv
import io
from pathlib import Path

from antismash_review.cli import main
from antismash_review.exporters.entity_tables import (
    DOMAIN_COLUMNS,
    GENE_COLUMNS,
    render_domain_tsv,
    render_gene_tsv,
)
from antismash_review.genbank import parse_genbank
from antismash_review.models import (
    Gene,
    Location,
    LocationPart,
    Record,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _minimal_record(**kwargs: object) -> Record:
    defaults: dict[str, object] = {
        "record_id": "REC1.1",
        "name": "REC1",
        "description": "test record",
        "length": 1000,
        "molecule_type": "DNA",
        "topology": "linear",
        "source_path": Path("/data/sample.gbk"),
        "source_sha256": "abc123sha",
        "antismash_version": "8.0.4",
        "organism": "synthetic",
        "taxonomy": [],
    }
    defaults.update(kwargs)
    return Record(**defaults)  # type: ignore[arg-type]


def _simple_location(
    start: int,
    end: int,
    strand: int = 1,
    partial: bool = False,
    cross_origin: bool = False,
) -> Location:
    part = LocationPart(start=start, end=end, strand=strand, fuzzy_start=partial)
    return Location(
        start=start,
        end=end,
        strand=strand,
        parts=(part,),
        cross_origin=cross_origin,
        original=f"{start}..{end}",
    )


def test_gene_tsv_exact_header_and_columns() -> None:
    record = parse_genbank(FIXTURES / "semantics.gb")[0]
    tsv = render_gene_tsv([record])
    reader = csv.reader(io.StringIO(tsv), delimiter="\t")
    rows = list(reader)
    assert rows[0] == list(GENE_COLUMNS)
    # semantics.gb has 3 CDS features (G1, G2, G3)
    assert len(rows) == 4
    for row in rows:
        assert len(row) == len(GENE_COLUMNS)


def test_domain_tsv_exact_header_and_columns() -> None:
    record = parse_genbank(FIXTURES / "semantics.gb")[0]
    tsv = render_domain_tsv([record])
    reader = csv.reader(io.StringIO(tsv), delimiter="\t")
    rows = list(reader)
    assert rows[0] == list(DOMAIN_COLUMNS)
    # semantics.gb has 2 aSDomain features (TIGR001, AMP-binding)
    assert len(rows) == 3
    for row in rows:
        assert len(row) == len(DOMAIN_COLUMNS)


def test_gene_tsv_values_zero_based_coords_booleans_and_lists() -> None:
    record = _minimal_record()
    record.genes.append(
        Gene(
            location=_simple_location(0, 30, strand=1, partial=True, cross_origin=False),
            locus_tag="LOCUS_1",
            gene="abcA",
            product="polyketide synthase\twith tab and\nnewline",
            protein_id="PROT1",
            translation=None,
            gene_kind="biosynthetic",
            gene_functions=[],
            ec_numbers=["2.3.1.-", "1.1.1.1"],
            db_xrefs=["GO:0009058"],
            notes=[],
            inference=[],
            region_numbers=[1, 2],
            candidate_cluster_numbers=[1],
            protocluster_numbers=[1],
            proto_core_numbers=[1],
            qualifiers={},
        )
    )
    tsv = render_gene_tsv([record])
    reader = csv.reader(io.StringIO(tsv), delimiter="\t")
    rows = list(reader)
    assert len(rows) == 2
    row = rows[1]
    assert row[0] == str(record.source_path)
    assert row[1] == "abc123sha"
    assert row[2] == "sample.gbk"
    assert row[3] == "REC1.1"
    assert row[4] == "LOCUS_1"
    assert row[5] == "abcA"
    assert row[6] == "polyketide synthase\twith tab and\nnewline"
    assert row[7] == "biosynthetic"
    assert row[8] == "0"
    assert row[9] == "30"
    assert row[10] == "1"
    assert row[11] == "true"
    assert row[12] == "false"
    assert row[13] == "2.3.1.-; 1.1.1.1"
    assert row[14] == "GO:0009058"
    assert row[15] == "1; 2"
    assert row[16] == "1"
    assert row[17] == "1"
    assert row[18] == "1"


def test_gene_tsv_empty_optional_fields() -> None:
    record = _minimal_record()
    record.genes.append(
        Gene(
            location=_simple_location(10, 50, strand=None, partial=False, cross_origin=False),
            locus_tag=None,
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
            region_numbers=[],
            candidate_cluster_numbers=[],
            protocluster_numbers=[],
            proto_core_numbers=[],
            qualifiers={},
        )
    )
    tsv = render_gene_tsv([record])
    reader = csv.reader(io.StringIO(tsv), delimiter="\t")
    rows = list(reader)
    row = rows[1]
    assert row[4] == ""  # locus_tag
    assert row[5] == ""  # gene
    assert row[6] == ""  # product
    assert row[7] == "unclassified"  # gene_kind
    assert row[8] == "10"  # start
    assert row[9] == "50"  # end
    assert row[10] == ""  # strand
    assert row[11] == "false"  # partial
    assert row[12] == "false"  # cross_origin
    assert row[13] == ""  # ec_numbers
    assert row[14] == ""  # db_xrefs
    assert row[15] == ""  # region_numbers
    assert row[16] == ""  # candidate_cluster_numbers
    assert row[17] == ""  # protocluster_numbers
    assert row[18] == ""  # proto_core_numbers


def test_domain_tsv_nrps_vs_non_nrps_and_optional_values() -> None:
    record = parse_genbank(FIXTURES / "semantics.gb")[0]
    tsv = render_domain_tsv([record])
    reader = csv.reader(io.StringIO(tsv), delimiter="\t")
    rows = list(reader)
    # Header + 2 domains
    assert len(rows) == 3

    # Row 1: TIGR001 (tool=tigrfam, not nrps_pks)
    tigr_row = rows[1]
    assert tigr_row[4] == "tigr_domain_1"  # domain_id
    assert tigr_row[5] == "TIGR001"  # name
    assert tigr_row[6] == "tigrfam"  # tool
    assert tigr_row[7] == "false"  # is_nrps_pks
    assert tigr_row[8] == ""  # locus_tag is not on this feature
    assert tigr_row[9] == "39"  # start (0-based of 40..50)
    assert tigr_row[10] == "50"  # end
    assert tigr_row[14] == ""  # protein_start None
    assert tigr_row[15] == ""  # protein_end None
    assert tigr_row[16] == ""  # score None
    assert tigr_row[17] == ""  # evalue None
    assert tigr_row[18] == "other"  # subtypes
    assert tigr_row[19] == ""  # specificity

    # Row 2: AMP-binding (tool=nrps_pks_domains, is_nrps_pks=true)
    nrps_row = rows[2]
    assert nrps_row[4] == "nrps_domain_2"
    assert nrps_row[5] == "AMP-binding"
    assert nrps_row[6] == "nrps_pks_domains"
    assert nrps_row[7] == "true"  # is_nrps_pks
    assert nrps_row[18] == "AMP-binding"  # subtypes
    assert nrps_row[19] == "Ser"  # specificity


def test_cli_gene_tsv_stdout_and_output_file(tmp_path: Path, capsys: object) -> None:
    assert main(["inspect", str(FIXTURES / "semantics.gb"), "--format", "gene-tsv"]) == 0
    captured = capsys.readouterr()  # type: ignore[union-attr]
    reader = csv.reader(io.StringIO(captured.out), delimiter="\t")
    rows = list(reader)
    assert rows[0] == list(GENE_COLUMNS)
    assert len(rows) == 4

    out_file = tmp_path / "genes.tsv"
    assert (
        main(
            [
                "inspect",
                str(FIXTURES / "semantics.gb"),
                "--format",
                "gene-tsv",
                "--output",
                str(out_file),
            ]
        )
        == 0
    )
    assert out_file.exists()
    reader_file = csv.reader(io.StringIO(out_file.read_text(encoding="utf-8")), delimiter="\t")
    assert list(reader_file) == rows


def test_cli_domain_tsv_stdout_and_output_file(tmp_path: Path, capsys: object) -> None:
    assert main(["inspect", str(FIXTURES / "semantics.gb"), "--format", "domain-tsv"]) == 0
    captured = capsys.readouterr()  # type: ignore[union-attr]
    reader = csv.reader(io.StringIO(captured.out), delimiter="\t")
    rows = list(reader)
    assert rows[0] == list(DOMAIN_COLUMNS)
    assert len(rows) == 3

    out_file = tmp_path / "domains.tsv"
    assert (
        main(
            [
                "inspect",
                str(FIXTURES / "semantics.gb"),
                "--format",
                "domain-tsv",
                "--output",
                str(out_file),
            ]
        )
        == 0
    )
    assert out_file.exists()
    reader_file = csv.reader(io.StringIO(out_file.read_text(encoding="utf-8")), delimiter="\t")
    assert list(reader_file) == rows


def test_cross_origin_gene_and_domain_tsv() -> None:
    record = parse_genbank(FIXTURES / "cross-origin.gb")[0]
    gene_tsv = render_gene_tsv([record])
    reader = csv.reader(io.StringIO(gene_tsv), delimiter="\t")
    rows = list(reader)
    assert len(rows) == 2
    row = rows[1]
    assert row[4] == "XORIGIN1"
    assert row[12] == "true"  # cross_origin
