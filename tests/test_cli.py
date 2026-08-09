from __future__ import annotations

from pathlib import Path

import pytest

from antismash_review.cli import load_review_records, main
from antismash_review.clusterblast import ClusterBlastParseError
from antismash_review.discovery import discover
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
        load_review_records(manifest, lenient=False)

    # lenient should add diagnostic
    records, _ = load_review_records(manifest, lenient=True)
    assert len(records) == 1
    codes = {d.code for d in records[0].diagnostics}
    assert "clusterblast_parse_failed" in codes
