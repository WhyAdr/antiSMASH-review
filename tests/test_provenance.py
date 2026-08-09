from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from antismash_review.compare import compare_records
from antismash_review.exporters.compare_json import dumps_comparison
from antismash_review.exporters.provenance import dumps_provenance, render_provenance_tsv
from antismash_review.genbank import parse_genbank
from antismash_review.models import AntiSmashProvenance
from tests.test_assemblyline import _record


def _write_record(
    path: Path, *, structured: dict[str, object] | None = None, comment: str | None = None
) -> Path:
    source = SeqRecord(Seq("ATGC"), id="PROV.1", name="PROV.1", description="provenance fixture")
    source.annotations["molecule_type"] = "DNA"
    if structured is not None:
        source.annotations["structured_comment"] = {"antiSMASH-Data": structured}
    if comment is not None:
        source.annotations["comment"] = comment
    SeqIO.write(source, path, "genbank")
    return path


def test_structured_comment_provenance_preserves_known_and_unknown_fields(tmp_path: Path) -> None:
    path = _write_record(
        tmp_path / "structured.gbk",
        structured={
            "Version": "8.0.4",
            "Run date": "2026-08-01 13:13:19",
            "Pfam version": "37.0",
            "Mystery": ["one", "two"],
        },
    )

    record = parse_genbank(path)[0]
    provenance = record.antismash_provenance

    assert record.antismash_version == "8.0.4"
    assert provenance.version == "8.0.4"
    assert provenance.run_date == "2026-08-01 13:13:19"
    assert provenance.pfam_version == "37.0"
    assert provenance.raw_fields["Mystery"] == ("one", "two")
    assert '"antismash_provenance"' not in dumps_provenance([record])
    assert '"Mystery": [' in dumps_provenance([record])
    assert '"antismash_provenance"' not in __import__("antismash_review").dumps_records([record])


def test_raw_comment_fallback_and_missing_values_remain_explicit(tmp_path: Path) -> None:
    comment = """##antiSMASH-Data-START##
Version:: 7.0.0
Unknown key:: first
Unknown key:: second
##antiSMASH-Data-END##"""
    record = parse_genbank(_write_record(tmp_path / "raw.gbk", comment=comment))[0]

    assert record.antismash_provenance.version == "7.0.0"
    assert record.antismash_provenance.raw_fields["Unknown key"] == ("first", "second")
    assert record.antismash_provenance.run_date is None
    assert record.antismash_provenance.pfam_version is None


def test_provenance_manifest_is_deduplicated_and_tsv_is_stable() -> None:
    left = _record([], [], products=["NRPS"])
    right = _record([], [], products=["PKS"])
    left.source_path = Path("run.gbk")
    right.source_path = Path("run.gbk")
    left.source_sha256 = right.source_sha256 = "sha"
    left.antismash_provenance = AntiSmashProvenance(version="8.0.4", run_date="today")
    right.antismash_provenance = AntiSmashProvenance(version="8.0.4", run_date="today")

    json_output = dumps_provenance([right, left])
    tsv_output = render_provenance_tsv([right, left])

    assert json_output.count('"source_path": "run.gbk"') == 1
    assert '"records": [' in json_output
    assert tsv_output.splitlines()[0].startswith("source_path\tsource_sha256")
    assert json_output == dumps_provenance([right, left])


def test_comparison_provenance_delta_is_tri_state_and_versioned() -> None:
    left = _record([], [], products=["NRPS"])
    right = _record([], [], products=["NRPS"])
    left.antismash_provenance = AntiSmashProvenance(
        version="8.0.3",
        pfam_version="36",
        raw_fields={"Version": ("8.0.3",), "Pfam version": ("36",)},
    )
    right.antismash_provenance = AntiSmashProvenance(
        version="8.0.4",
        pfam_version="36",
        raw_fields={"Version": ("8.0.4",), "Pfam version": ("36",)},
    )
    result = compare_records(
        [left],
        [right],
        left_input=Path("left"),
        right_input=Path("right"),
    )
    delta = result.matched[0].provenance_delta

    assert delta is not None
    assert delta.antismash_version_changed is True
    assert delta.pfam_version_changed is False
    assert delta.detection_rule_set_changed is None
    assert delta.differing_raw_fields == ("Version",)
    assert '"schema_version": "0.2.0"' in dumps_comparison(result)

    unknown_left = _record([], [], products=["NRPS"])
    unknown_right = _record([], [], products=["NRPS"])
    unknown = compare_records(
        [unknown_left],
        [unknown_right],
        left_input=Path("left"),
        right_input=Path("right"),
    )
    unknown_delta = unknown.matched[0].provenance_delta
    assert unknown_delta is not None
    assert unknown_delta.antismash_version_changed is None
