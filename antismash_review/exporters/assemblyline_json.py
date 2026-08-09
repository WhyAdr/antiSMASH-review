"""JSON export for evidence-only assembly-line predictions."""

from __future__ import annotations

import json
from dataclasses import asdict
from typing import Any

from antismash_review._version import __version__
from antismash_review.assemblyline import predict_assembly_lines
from antismash_review.models import Record
from antismash_review.schema import ASSEMBLYLINE_SCHEMA_NAME, ASSEMBLYLINE_SCHEMA_VERSION


def dumps_assembly_lines(records: list[Record]) -> str:
    document: dict[str, Any] = {
        "schema_name": ASSEMBLYLINE_SCHEMA_NAME,
        "schema_version": ASSEMBLYLINE_SCHEMA_VERSION,
        "parser_version": __version__,
        "predictions": [
            asdict(prediction)
            for record in records
            for prediction in predict_assembly_lines(record)
        ],
    }
    return json.dumps(document, indent=2, sort_keys=True) + "\n"
