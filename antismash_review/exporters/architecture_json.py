"""JSON export for scoped domain-architecture assessments."""

from __future__ import annotations

import json
from dataclasses import asdict

from antismash_review._version import __version__
from antismash_review.architecture import assess_architecture
from antismash_review.models import Record


def dumps_architecture(records: list[Record]) -> str:
    document = {
        "schema_name": "antismash-review-architecture",
        "schema_version": "0.1.0",
        "parser_version": __version__,
        "assessments": [
            asdict(assessment) for record in records for assessment in assess_architecture(record)
        ],
    }
    return json.dumps(document, indent=2, sort_keys=True) + "\n"
