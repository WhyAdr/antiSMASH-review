from __future__ import annotations

import json
from dataclasses import asdict
from enum import Enum
from pathlib import Path
from typing import Any

from antismash_review._version import __version__
from antismash_review.compare import ComparisonResult
from antismash_review.schema import COMPARISON_SCHEMA_NAME, COMPARISON_SCHEMA_VERSION


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def dumps_comparison(result: ComparisonResult) -> str:
    document = {
        "schema_name": COMPARISON_SCHEMA_NAME,
        "schema_version": COMPARISON_SCHEMA_VERSION,
        "parser_version": __version__,
        "comparison": asdict(result),
    }
    return json.dumps(document, indent=2, sort_keys=True, default=_json_default) + "\n"
