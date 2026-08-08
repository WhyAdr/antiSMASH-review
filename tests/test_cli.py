import json
import shutil
from pathlib import Path

import pytest

from antismash_review.cli import main

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "semantics.gb"


def test_cli_json_is_versioned_and_serializable(capsys: object) -> None:
    assert main(["inspect", str(FIXTURE), "--format", "json"]) == 0
    output = capsys.readouterr().out  # type: ignore[union-attr]
    document = json.loads(output)
    assert document["schema_name"] == "antismash-review"
    assert document["records"][0]["antismash_version"] == "8.0.4"
    assert document["records"][0]["source_path"]


def test_cli_tsv_has_compact_record_columns(capsys: object) -> None:
    assert main(["inspect", str(FIXTURE), "--format", "tsv"]) == 0
    output = capsys.readouterr().out  # type: ignore[union-attr]
    lines = output.splitlines()
    assert lines[0].startswith("filename\trecord_id\tregion_products")
    assert "\tall_domains\tmodules\traw_pfam_hits\tdeduplicated_pfam_hits\t" in lines[0]
    assert "hybrid-a; hybrid-b" in lines[1]


def test_cli_can_write_markdown(tmp_path: Path) -> None:
    output = tmp_path / "review.md"
    assert main(["inspect", str(FIXTURE), "--output", str(output)]) == 0
    assert "# antiSMASH review" in output.read_text(encoding="utf-8")


def test_cli_refuses_to_overwrite_input(tmp_path: Path, capsys: object) -> None:
    source = tmp_path / "source.gb"
    shutil.copyfile(FIXTURE, source)
    original = source.read_bytes()

    assert main(["inspect", str(source), "--output", str(source)]) == 2
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "refusing to overwrite input file" in captured.err
    assert source.read_bytes() == original


def test_cli_reports_output_write_failure(tmp_path: Path, capsys: object) -> None:
    output = tmp_path / "missing" / "review.md"
    assert main(["inspect", str(FIXTURE), "--output", str(output)]) == 2
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "could not write" in captured.err


def test_cli_reports_native_json_as_unsupported(tmp_path: Path, capsys: object) -> None:
    source = tmp_path / "results.json"
    source.write_text("{}", encoding="utf-8")
    assert main(["inspect", str(source)]) == 2
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "native antiSMASH JSON cannot yet provide the review record model" in captured.err


def test_cli_no_command_shows_help() -> None:
    """Calling the CLI with no arguments should fail (required subcommand)."""
    with pytest.raises(SystemExit) as exc_info:
        main([])
    assert exc_info.value.code == 2


def test_cli_discovery_error_for_nonexistent_path(capsys: object) -> None:
    """CLI returns status 2 and stderr message for nonexistent input path."""
    result = main(["inspect", "/nonexistent/path/to/nothing.gbk"])
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert result == 2
    assert "error:" in captured.err.casefold()


def test_cli_parse_error_returns_status_2(tmp_path: Path, capsys: object) -> None:
    """CLI returns status 2 when GenBank parsing fails on corrupt content."""
    corrupt = tmp_path / "corrupt.gb"
    corrupt.write_text("THIS IS NOT GENBANK FORMAT", encoding="utf-8")
    result = main(["inspect", str(corrupt)])
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert result == 2
    assert "error:" in captured.err.casefold()


def test_cli_json_only_directory_rejects(tmp_path: Path, capsys: object) -> None:
    """A directory containing only JSON files returns status 2."""
    (tmp_path / "result.json").write_text("{}", encoding="utf-8")
    result = main(["inspect", str(tmp_path)])
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert result == 2
    assert "native antiSMASH JSON" in captured.err or "no GenBank" in captured.err
