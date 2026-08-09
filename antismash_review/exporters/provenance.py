"""Source-level provenance manifest exports."""

from __future__ import annotations

import csv
import io
import json
from dataclasses import asdict
from pathlib import Path

from antismash_review._version import __version__
from antismash_review.models import AntiSmashProvenance, Record
from antismash_review.schema import PROVENANCE_SCHEMA_NAME, PROVENANCE_SCHEMA_VERSION


def _merge_provenance(values: list[AntiSmashProvenance]) -> AntiSmashProvenance:
    def same(attribute: str) -> str | None:
        candidates = {getattr(value, attribute) for value in values}
        return next(iter(candidates)) if len(candidates) == 1 else None

    raw_keys = {key for value in values for key in value.raw_fields}
    raw_fields = {
        key: tuple(
            dict.fromkeys(item for value in values for item in value.raw_fields.get(key, ()))
        )
        for key in sorted(raw_keys)
    }
    database_keys = {key for value in values for key in value.database_versions}
    database_versions = {
        key: values[0].database_versions[key]
        for key in sorted(database_keys)
        if len({value.database_versions.get(key) for value in values}) == 1
        and key in values[0].database_versions
    }
    return AntiSmashProvenance(
        version=same("version"),
        run_date=same("run_date"),
        pfam_version=same("pfam_version"),
        detection_rule_set_version=same("detection_rule_set_version"),
        database_versions=database_versions,
        raw_fields=raw_fields,
    )


def _manifest_inputs(records: list[Record]) -> list[dict[str, object]]:
    grouped: dict[tuple[Path, str], list[Record]] = {}
    for record in records:
        grouped.setdefault((record.source_path, record.source_sha256), []).append(record)

    inputs: list[dict[str, object]] = []
    for (source_path, source_sha256), source_records in sorted(
        grouped.items(), key=lambda item: (str(item[0][0]), item[0][1])
    ):
        provenance = _merge_provenance([record.antismash_provenance for record in source_records])
        record_ids = sorted({record.record_id for record in source_records})
        inputs.append(
            {
                "source_path": str(source_path),
                "source_sha256": source_sha256,
                "records": record_ids,
                "antismash": asdict(provenance),
            }
        )
    return inputs


def dumps_provenance(records: list[Record]) -> str:
    document = {
        "schema_name": PROVENANCE_SCHEMA_NAME,
        "schema_version": PROVENANCE_SCHEMA_VERSION,
        "review_tool": {"name": "antismash-review", "version": __version__},
        "inputs": _manifest_inputs(records),
    }
    return json.dumps(document, indent=2, sort_keys=True) + "\n"


PROVENANCE_COLUMNS = (
    "source_path",
    "source_sha256",
    "records",
    "antismash_version",
    "run_date",
    "pfam_version",
    "detection_rule_set_version",
    "database_versions",
    "raw_fields",
)


def render_provenance_tsv(records: list[Record]) -> str:
    document = json.loads(dumps_provenance(records))
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    writer.writerow(PROVENANCE_COLUMNS)
    for item in document["inputs"]:
        antismash = item["antismash"]
        writer.writerow(
            (
                item["source_path"],
                item["source_sha256"],
                ";".join(item["records"]),
                antismash["version"] or "",
                antismash["run_date"] or "",
                antismash["pfam_version"] or "",
                antismash["detection_rule_set_version"] or "",
                json.dumps(antismash["database_versions"], sort_keys=True),
                json.dumps(antismash["raw_fields"], sort_keys=True),
            )
        )
    return output.getvalue()
