from pathlib import Path

import pytest

from antismash_review.discovery import discover


def test_discovery_separates_aggregate_and_region_inputs(tmp_path: Path) -> None:
    aggregate = tmp_path / "final.gbk"
    region = tmp_path / "sample.region001.gbk"
    ignored = tmp_path / "notes.txt"
    for path in (aggregate, region, ignored):
        path.write_text("", encoding="utf-8")

    manifest = discover(tmp_path)
    assert manifest.aggregate_genbanks == (aggregate.resolve(),)
    assert manifest.region_genbanks == (region.resolve(),)
    assert manifest.ignored_files == (ignored.resolve(),)


def test_recursive_discovery_is_deterministic(tmp_path: Path) -> None:
    nested = tmp_path / "run" / "nested"
    nested.mkdir(parents=True)
    path = nested / "sample.gbff"
    path.write_text("", encoding="utf-8")
    manifest = discover(tmp_path, recursive=True)
    assert manifest.aggregate_genbanks == (path.resolve(),)


def test_region_classification_is_case_insensitive(tmp_path: Path) -> None:
    path = tmp_path / "sample.REGION001.GBK"
    path.write_text("", encoding="utf-8")
    manifest = discover(path)
    assert manifest.region_genbanks == (path.resolve(),)
    assert not manifest.aggregate_genbanks


def test_unsupported_single_file_raises_error(tmp_path: Path) -> None:
    """A single file with an unsupported extension should raise ValueError."""
    unsupported = tmp_path / "data.csv"
    unsupported.write_text("a,b,c", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported input"):
        discover(unsupported)


def test_nonexistent_directory_raises_error() -> None:
    """A nonexistent directory should raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        discover(Path("/nonexistent/antismash/output"))
