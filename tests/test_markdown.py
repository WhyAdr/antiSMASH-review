from __future__ import annotations

from pathlib import Path

from antismash_review.exporters.markdown import _escape_cell, _text, render_records
from antismash_review.genbank import parse_genbank
from antismash_review.models import (
    ClusterBlastHit,
    ClusterBlastPairing,
    ClusterBlastResult,
    Diagnostic,
    Severity,
)
from tests.fixtures.build_fixture import write_synthetic_genbank


def test_escape_cell_and_text_helpers() -> None:
    assert _text(None) == "not reported"
    assert _text("") == "not reported"
    assert _text("8.0.0") == "8.0.0"

    assert _escape_cell(None) == ""
    assert _escape_cell("abc|def") == r"abc\|def"
    assert _escape_cell("abc\r\ndef\tghi") == "abc def ghi"


def test_render_records_empty() -> None:
    rendered = render_records([])
    assert rendered == "# antiSMASH review\n"


def test_render_records_with_input_diagnostics(tmp_path: Path) -> None:
    gbk = write_synthetic_genbank(tmp_path / "synthetic.gbk")
    records = parse_genbank(gbk)

    diag = Diagnostic(
        code="clusterblast_parse_failed",
        severity=Severity.WARNING,
        message="Sidecar file was corrupted",
        source="some/path.txt",
    )

    rendered = render_records(records, input_diagnostics=[diag])
    assert "## Input diagnostics" in rendered
    assert "- `warning` `clusterblast_parse_failed`: Sidecar file was corrupted" in rendered
    assert f"## `{records[0].record_id}`" in rendered


def test_render_records_synthetic_genbank(tmp_path: Path) -> None:
    gbk = write_synthetic_genbank(tmp_path / "synthetic.gbk")
    records = parse_genbank(gbk)

    rendered1 = render_records(records)
    rendered2 = render_records(records)
    assert rendered1 == rendered2

    assert "# antiSMASH review" in rendered1
    assert "## `SYNTH.1`" in rendered1
    assert "### Products" in rendered1
    assert "### Diagnostics" in rendered1


def test_render_records_clusterblast_table_and_pairings(tmp_path: Path) -> None:
    gbk = write_synthetic_genbank(tmp_path / "synthetic.gbk")
    records = parse_genbank(gbk)

    pairing = ClusterBlastPairing(
        query_gene="GENE_1",
        subject_gene="SUBJ_1",
        percent_identity=90.5,
        blast_score=150.0,
        percent_coverage=85.0,
        evalue=1e-25,
        subject_protein_id="PROT_123",
    )
    hit = ClusterBlastHit(
        rank=1,
        accession="BGC0001234",
        description="Test BGC with | bar",
        cluster_type="NRPS",
        num_hits=1,
        core_gene_hits=1,
        blast_score=150.0,
        synteny_score=2,
        core_bonus=1,
        similarity=85,
        pairings=[pairing],
    )
    cb_result = ClusterBlastResult(
        record_id="SYNTH.1",
        region_number=1,
        search_type="clusterblast",
        total_hits=1,
        rankings=[hit],
        source_path=Path("test/clusterblast.txt"),
        source_sha256="abc123sha",
        source_format="txt",
    )
    records[0].clusterblast_results.append(cb_result)

    rendered = render_records(records)
    assert "### ClusterBlast" in rendered
    assert "#### Region 1 — ClusterBlast" in rendered
    assert "BGC0001234" in rendered
    assert r"Test BGC with \| bar" in rendered
    assert "85" in rendered
