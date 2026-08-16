from __future__ import annotations

import json
from dataclasses import asdict
from enum import Enum
from pathlib import Path
from typing import Any

from antismash_review._version import __version__
from antismash_review.models import Record
from antismash_review.review import review_record
from antismash_review.schema import RECORD_SCHEMA_NAME, RECORD_SCHEMA_VERSION


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def dumps_records(records: list[Record]) -> str:
    serialized_records = []
    for record in records:
        serialized = asdict(record)
        # Schema 0.3.0: asdict() now deliberately includes antismash_provenance.
        # This exposes version, run_date, pfam_version, database_versions, raw_fields.
        serialized_records.append(serialized)
    document = {
        "schema_name": RECORD_SCHEMA_NAME,
        "schema_version": RECORD_SCHEMA_VERSION,
        "parser_version": __version__,
        "records": serialized_records,
        "diagnostics": [
            asdict(diagnostic) for record in records for diagnostic in review_record(record)
        ],
    }
    return json.dumps(document, indent=2, sort_keys=True, default=_json_default) + "\n"
