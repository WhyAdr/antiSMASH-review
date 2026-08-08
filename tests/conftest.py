from __future__ import annotations

import os
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture
def sm_zpg19_dir() -> Path:
    configured = os.environ.get("ANTISMASH_SM_ZPG19_DIR")
    candidates = [
        Path(configured).expanduser() if configured else None,
        REPOSITORY_ROOT / "SM-ZPG19--NOACCESSION-antismash",
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_dir():
            return candidate.resolve()
    pytest.skip("private SM-ZPG19 antiSMASH result directory is unavailable")
