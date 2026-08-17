from __future__ import annotations

from pathlib import Path

import pytest

from antismash_review.cli import load_review_records, main
from antismash_review.clusterblast import (
    ClusterBlastParseError,
    attach_clusterblast_results,
    merge_clusterblast_results,
    parse_clusterblast_text,
)
from antismash_review.discovery import discover
from antismash_review.genbank import parse_genbank
from antismash_review.loading import load_review_input
from antismash_review.models import Diagnostic
from tests.fixtures.build_fixture import write_synthetic_genbank


@pytest.mark.parametrize(
    "output_format",
    ["markdown", "json", "tsv", "gene-tsv", "domain-tsv", "clusterblast-tsv"],
)
def test_cli_inspect_dispatches_every_public_format(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    output_format: str,
) -> None:
    input_path = write_synthetic_genbank(tmp_path / "synthetic.gbk")

    assert main(["inspect", str(input_path), "--format", output_format]) == 0
    output = capsys.readouterr().out
    assert output
    if output_format == "json":
        assert '"records"' in output
    elif output_format == "markdown":
        assert "# antiSMASH review" in output
    else:
        assert output.splitlines()[0]


@pytest.mark.parametrize(
    "output_format",
    [
        "assemblyline-tsv",
        "assemblyline-json",
        "assemblyline-markdown",
        "architecture-json",
        "gff3",
        "bed",
        "provenance-json",
        "provenance-tsv",
    ],
)
def test_cli_inspect_dispatches_assemblyline_formats(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    output_format: str,
) -> None:
    input_path = write_synthetic_genbank(tmp_path / "synthetic.gbk")

    assert main(["inspect", str(input_path), "--format", output_format]) == 0
    output = capsys.readouterr().out
    assert output


def test_cli_inspect_writes_output_file(
    tmp_path: Path,
) -> None:
    input_path = write_synthetic_genbank(tmp_path / "synthetic.gbk")
    output_path = tmp_path / "output.json"

    assert main(["inspect", str(input_path), "--format", "json", "--output", str(output_path)]) == 0
    assert output_path.is_file()
    assert '"records"' in output_path.read_text(encoding="utf-8")


def test_cli_inspect_refuses_to_overwrite_input(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    input_path = write_synthetic_genbank(tmp_path / "synthetic.gbk")

    assert main(["inspect", str(input_path), "--output", str(input_path)]) == 2
    assert "refusing to overwrite input file" in capsys.readouterr().err


def test_cli_inspect_nonexistent_input_returns_error(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert main(["inspect", str(tmp_path / "nonexistent.gbk")]) == 2
    assert "error:" in capsys.readouterr().err


def test_cli_coordinate_mode_requires_explicit_shared_coordinates() -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["compare", "left.gbk", "right.gbk", "--match-by", "coordinate_overlap"])
    assert exc_info.value.code == 2


def test_cli_coordinate_mode_validates_overlap_range() -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "compare",
                "left.gbk",
                "right.gbk",
                "--match-by",
                "coordinate_overlap",
                "--assume-shared-coordinate-system",
                "--min-reciprocal-overlap",
                "1.5",
            ]
        )
    assert exc_info.value.code == 2


def test_cli_shared_coordinates_flag_rejected_without_coordinate_overlap() -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "compare",
                "left.gbk",
                "right.gbk",
                "--match-by",
                "record_id",
                "--assume-shared-coordinate-system",
            ]
        )
    assert exc_info.value.code == 2


def test_cli_compare_runs_markdown_and_json(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    left_path = write_synthetic_genbank(tmp_path / "left.gbk")
    right_path = write_synthetic_genbank(tmp_path / "right.gbk")

    # Markdown format
    assert main(["compare", str(left_path), str(right_path), "--format", "markdown"]) == 0
    out = capsys.readouterr().out
    assert "Comparison" in out or "SYNTH.1" in out

    # JSON format to file
    out_json = tmp_path / "comparison.json"
    assert (
        main(
            [
                "compare",
                str(left_path),
                str(right_path),
                "--format",
                "json",
                "--output",
                str(out_json),
            ]
        )
        == 0
    )
    assert out_json.is_file()
    assert "antismash-review-comparison" in out_json.read_text(encoding="utf-8")


def test_cli_compare_refuses_to_overwrite_input(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    left_path = write_synthetic_genbank(tmp_path / "left.gbk")
    right_path = write_synthetic_genbank(tmp_path / "right.gbk")

    assert main(["compare", str(left_path), str(right_path), "--output", str(left_path)]) == 2
    assert "refusing to overwrite input file" in capsys.readouterr().err


def test_cli_compare_error_on_missing_input(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    left_path = write_synthetic_genbank(tmp_path / "left.gbk")
    assert main(["compare", str(left_path), str(tmp_path / "missing.gbk")]) == 2
    assert "error:" in capsys.readouterr().err


def test_cli_compare_coordinate_overlap_dispatch(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    left_path = write_synthetic_genbank(tmp_path / "left.gbk")
    right_path = write_synthetic_genbank(tmp_path / "right.gbk")

    assert (
        main(
            [
                "compare",
                str(left_path),
                str(right_path),
                "--match-by",
                "coordinate_overlap",
                "--assume-shared-coordinate-system",
                "--min-reciprocal-overlap",
                "0.80",
            ]
        )
        == 0
    )
    out = capsys.readouterr().out
    assert "SYNTH.1" in out


def test_load_review_records_json_only_error(tmp_path: Path) -> None:
    json_file = tmp_path / "results.json"
    json_file.write_text('{"records": []}', encoding="utf-8")
    manifest = discover(tmp_path)
    with pytest.raises(ValueError, match="native antiSMASH JSON cannot yet provide"):
        load_review_records(manifest, lenient=False)


def test_load_review_records_empty_manifest_error(tmp_path: Path) -> None:
    manifest = discover(tmp_path)
    with pytest.raises(ValueError, match="no GenBank input found"):
        load_review_records(manifest, lenient=False)


def test_load_review_records_lenient_sidecar_diagnostic(tmp_path: Path) -> None:
    write_synthetic_genbank(tmp_path / "synthetic.gbk")
    cb_dir = tmp_path / "clusterblast"
    cb_dir.mkdir()
    bad_cb = cb_dir / "contig_1_c1.txt"
    bad_cb.write_text("invalid clusterblast content", encoding="utf-8")

    manifest = discover(tmp_path)
    # strict should raise
    with pytest.raises(ClusterBlastParseError):
        load_review_input(manifest, lenient=False)

    # lenient should add diagnostic to input_diagnostics
    loaded = load_review_input(manifest, lenient=True)
    assert len(loaded.records) == 1
    assert any(d.code == "clusterblast_parse_failed" for d in loaded.input_diagnostics)
    assert not any(d.code == "clusterblast_parse_failed" for d in loaded.records[0].diagnostics)


def test_lenient_loading_retains_valid_sidecars_alongside_corrupt_ones(
    tmp_path: Path,
) -> None:
    write_synthetic_genbank(tmp_path / "synthetic.gbk")
    cb_dir = tmp_path / "clusterblast"
    cb_dir.mkdir()
    valid_cb = cb_dir / "contig_1_c1.txt"
    valid_cb.write_text(
        "ClusterBlast scores for SYNTH.1\n"
        "Significant hits:\n"
        "1. SYNTH-HIT-1\tSynthetic hit\n"
        "Details:\n"
        "1. SYNTH-HIT-1\n"
        "Type: NRPS\n"
        "Number of proteins with BLAST hits to this cluster: 1\n"
        "Cumulative BLAST score: 12.5\n"
        "Table of Blast hits\n"
        "SYN_CDS_1\tSYNTH_SUBJECT\t90.0\t12.5\t80.0\t1e-10\n"
        ">>\n",
        encoding="utf-8",
    )

    kcb_dir = tmp_path / "knownclusterblast"
    kcb_dir.mkdir()
    bad_kcb = kcb_dir / "contig_1_c1.txt"
    bad_kcb.write_text("corrupt knownclusterblast content\n", encoding="utf-8")

    manifest = discover(tmp_path)

    # strict mode must raise
    with pytest.raises(ClusterBlastParseError):
        load_review_input(manifest, lenient=False)

    # lenient mode must retain the valid clusterblast result and emit diagnostic for the corrupt one
    loaded = load_review_input(manifest, lenient=True)
    assert len(loaded.records) == 1
    assert len(loaded.records[0].clusterblast_results) == 1
    assert loaded.records[0].clusterblast_results[0].search_type == "clusterblast"
    assert loaded.records[0].clusterblast_results[0].rankings[0].accession == "SYNTH-HIT-1"

    assert any(d.code == "clusterblast_parse_failed" for d in loaded.input_diagnostics)
    assert not any(d.code == "clusterblast_parse_failed" for d in loaded.records[0].diagnostics)


def test_lenient_loading_retains_valid_sidecars_when_one_target_is_unattachable(
    tmp_path: Path,
) -> None:
    write_synthetic_genbank(tmp_path / "synthetic.gbk")
    cb_dir = tmp_path / "clusterblast"
    cb_dir.mkdir()
    # valid target for SYNTH.1 region 1
    valid_cb = cb_dir / "contig_1_c1.txt"
    valid_cb.write_text(
        "ClusterBlast scores for SYNTH.1\n"
        "Significant hits:\n"
        "1. SYNTH-HIT-1\tSynthetic hit\n"
        "Details:\n"
        "1. SYNTH-HIT-1\n"
        "Type: NRPS\n"
        "Number of proteins with BLAST hits to this cluster: 1\n"
        "Cumulative BLAST score: 12.5\n"
        "Table of Blast hits\n"
        "SYN_CDS_1\tSYNTH_SUBJECT\t90.0\t12.5\t80.0\t1e-10\n"
        ">>\n",
        encoding="utf-8",
    )
    # unattachable target for SYNTH.1 region 99 (no region 99 exists)
    unattachable_cb = cb_dir / "contig_1_c99.txt"
    unattachable_cb.write_text(
        "ClusterBlast scores for SYNTH.1\n"
        "Significant hits:\n"
        "1. SYNTH-HIT-99\tSynthetic hit 99\n"
        "Details:\n"
        "1. SYNTH-HIT-99\n"
        "Type: NRPS\n"
        "Number of proteins with BLAST hits to this cluster: 1\n"
        "Cumulative BLAST score: 12.5\n"
        "Table of Blast hits\n"
        "SYN_CDS_1\tSYNTH_SUBJECT\t90.0\t12.5\t80.0\t1e-10\n"
        ">>\n",
        encoding="utf-8",
    )

    manifest = discover(tmp_path)

    # strict mode must fail before mutating records
    with pytest.raises(ClusterBlastParseError, match="expected one GenBank target"):
        load_review_records(manifest, lenient=False)

    # lenient mode attaches the valid result and emits clusterblast_attach_failed for the other
    records, _ = load_review_records(manifest, lenient=True)
    assert len(records) == 1
    assert len(records[0].clusterblast_results) == 1
    assert records[0].clusterblast_results[0].region_number == 1
    assert records[0].clusterblast_results[0].rankings[0].accession == "SYNTH-HIT-1"

    codes = {d.code for d in records[0].diagnostics}
    assert "clusterblast_attach_failed" in codes


def test_lenient_loading_retains_valid_sidecars_when_duplicate_results_occur(
    tmp_path: Path,
) -> None:
    # Direct testing of merge_clusterblast_results and attach_clusterblast_results in lenient mode
    records = parse_genbank(write_synthetic_genbank(tmp_path / "synthetic.gbk"))
    cb_file = tmp_path / "contig_1_c1.txt"
    cb_file.write_text(
        "ClusterBlast scores for SYNTH.1\n"
        "Significant hits:\n"
        "1. SYNTH-HIT-1\tSynthetic hit\n"
        "Details:\n"
        "1. SYNTH-HIT-1\n"
        "Type: NRPS\n"
        "Number of proteins with BLAST hits to this cluster: 1\n"
        "Cumulative BLAST score: 12.5\n"
        "Table of Blast hits\n"
        "SYN_CDS_1\tSYNTH_SUBJECT\t90.0\t12.5\t80.0\t1e-10\n"
        ">>\n",
        encoding="utf-8",
    )
    res1 = parse_clusterblast_text(cb_file, search_type="clusterblast")
    res2 = parse_clusterblast_text(cb_file, search_type="clusterblast")

    # Duplicate text results
    diagnostics: list[Diagnostic] = []
    merged = merge_clusterblast_results([res1, res2], [], lenient=True, diagnostics=diagnostics)
    assert len(merged) == 1
    assert len(diagnostics) == 1
    assert diagnostics[0].code == "clusterblast_duplicate_result"

    attach_clusterblast_results(records, merged, lenient=True)
    assert len(records[0].clusterblast_results) == 1
    assert records[0].clusterblast_results[0].rankings[0].accession == "SYNTH-HIT-1"


def test_lenient_loading_diagnostic_provenance_multi_record(tmp_path: Path) -> None:
    rec1_path = write_synthetic_genbank(tmp_path / "rec1.gbk")
    rec1 = parse_genbank(rec1_path)[0]
    rec2_path = write_synthetic_genbank(tmp_path / "rec2.gbk")
    rec2 = parse_genbank(rec2_path)[0]
    rec2.record_id = "RECORD_B"
    rec2.name = "RECORD_B"
    records = [rec1, rec2]

    cb_file = tmp_path / "contig_1_c99.txt"
    cb_file.write_text(
        "ClusterBlast scores for RECORD_B\n"
        "Significant hits:\n"
        "1. HIT-B\tHit description\n"
        "Details:\n"
        "1. HIT-B\n"
        "Type: NRPS\n"
        "Number of proteins with BLAST hits to this cluster: 1\n"
        "Cumulative BLAST score: 10.0\n"
        "Table of Blast hits\n"
        "SYN_CDS_1\tSUBJ_1\t90.0\t10.0\t80.0\t1e-10\n"
        ">>\n",
        encoding="utf-8",
    )
    result = parse_clusterblast_text(cb_file, search_type="clusterblast")

    attach_clusterblast_results(records, [result], lenient=True)

    # Diagnostic must be attached to RECORD_B (rec2), not rec1
    assert not any(d.code == "clusterblast_attach_failed" for d in rec1.diagnostics)
    attach_diags = [d for d in rec2.diagnostics if d.code == "clusterblast_attach_failed"]
    assert len(attach_diags) == 1
    assert attach_diags[0].record_id == "RECORD_B"


def test_lenient_loading_diagnostic_isolation_multi_record(tmp_path: Path) -> None:
    write_synthetic_genbank(tmp_path / "rec1.gbk")
    # second record in separate file
    write_synthetic_genbank(tmp_path / "rec2.gbk")
    # Update rec2 gbk file to have RECORD_B as locus/accession
    rec2_content = (tmp_path / "rec2.gbk").read_text(encoding="utf-8")
    rec2_content = rec2_content.replace("SYNTH.1", "RECORD_B")
    (tmp_path / "rec2.gbk").write_text(rec2_content, encoding="utf-8")

    # unattachable sidecar for RECORD_B region 99
    cb_dir = tmp_path / "clusterblast"
    cb_dir.mkdir()
    cb_file = cb_dir / "contig_1_c99.txt"
    cb_file.write_text(
        "ClusterBlast scores for RECORD_B\n"
        "Significant hits:\n"
        "1. HIT-B\tHit description\n"
        "Details:\n"
        "1. HIT-B\n"
        "Type: NRPS\n"
        "Number of proteins with BLAST hits to this cluster: 1\n"
        "Cumulative BLAST score: 10.0\n"
        "Table of Blast hits\n"
        "SYN_CDS_1\tSUBJ_1\t90.0\t10.0\t80.0\t1e-10\n"
        ">>\n",
        encoding="utf-8",
    )

    manifest = discover(tmp_path)
    loaded = load_review_input(manifest, lenient=True)

    rec_a = next(r for r in loaded.records if r.record_id == "SYNTH.1")
    rec_b = next(r for r in loaded.records if r.record_id == "RECORD_B")

    # rec_a must NOT have any attach failed diagnostics
    assert not any(d.code == "clusterblast_attach_failed" for d in rec_a.diagnostics)
    # rec_b must have the attach failed diagnostic
    assert any(d.code == "clusterblast_attach_failed" for d in rec_b.diagnostics)
