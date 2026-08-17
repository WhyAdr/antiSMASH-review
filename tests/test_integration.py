from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from antismash_review.discovery import discover
from antismash_review.loading import load_review_input


def test_integration_manifest_harness() -> None:
    manifest_path_env = os.environ.get("ANTISMASH_REVIEW_INTEGRATION_MANIFEST")
    if not manifest_path_env:
        pytest.skip(
            "ANTISMASH_REVIEW_INTEGRATION_MANIFEST not set; skipping real-data integration suite."
        )

    manifest_path = Path(manifest_path_env)
    if not manifest_path.is_file():
        pytest.skip(f"Integration manifest file does not exist: {manifest_path}")

    content = manifest_path.read_text(encoding="utf-8")
    data = json.loads(content)

    for entry in data.get("cases", []):
        input_dir = Path(entry["path"])
        if not input_dir.exists():
            continue
        discovered = discover(input_dir, recursive=entry.get("recursive", False))
        loaded = load_review_input(discovered, lenient=entry.get("lenient", False))

        if "expected_records_count" in entry:
            assert len(loaded.records) == entry["expected_records_count"]
        if "expected_record_ids" in entry:
            assert [r.record_id for r in loaded.records] == entry["expected_record_ids"]
