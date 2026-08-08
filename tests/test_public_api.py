from __future__ import annotations

from importlib import resources
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import antismash_review
from antismash_review import (
    Diagnostic,
    GenBankParseError,
    Record,
    Severity,
    dumps_records,
    parse_genbank,
    render_clusterblast_tsv,
    render_domain_tsv,
    render_gene_tsv,
    render_records,
    render_tsv,
    review_record,
)

EXPECTED_PUBLIC_NAMES = {
    "Diagnostic",
    "GenBankParseError",
    "Record",
    "Severity",
    "__version__",
    "dumps_records",
    "parse_genbank",
    "render_clusterblast_tsv",
    "render_domain_tsv",
    "render_gene_tsv",
    "render_records",
    "render_tsv",
    "review_record",
}


def _write_minimal_genbank(tmp_path: Path) -> Path:
    source = SeqRecord(Seq("ATGC"), id="api_record", name="api_record", description="API test")
    source.annotations["molecule_type"] = "DNA"
    path = tmp_path / "api.gb"
    SeqIO.write(source, path, "genbank")
    return path


def test_top_level_public_api_exports_expected_names() -> None:
    assert set(antismash_review.__all__) == EXPECTED_PUBLIC_NAMES
    assert antismash_review.parse_genbank is parse_genbank
    assert antismash_review.GenBankParseError is GenBankParseError
    assert antismash_review.Record is Record
    assert antismash_review.Diagnostic is Diagnostic
    assert antismash_review.Severity is Severity


def test_top_level_parse_review_export_workflow(tmp_path: Path) -> None:
    records = parse_genbank(_write_minimal_genbank(tmp_path))

    assert len(records) == 1
    assert isinstance(records[0], Record)
    assert records[0].record_id == "api_record"
    assert review_record(records[0]) == []
    assert '"record_id": "api_record"' in dumps_records(records)
    assert "# antiSMASH review" in render_records(records)
    assert "api_record" in render_tsv(records)
    assert render_gene_tsv(records).startswith("source_path\t")
    assert render_domain_tsv(records).startswith("source_path\t")
    assert render_clusterblast_tsv(records).startswith("record_id\t")


def test_pep561_marker_is_available_as_package_data() -> None:
    marker = resources.files("antismash_review").joinpath("py.typed")
    assert marker.is_file()
