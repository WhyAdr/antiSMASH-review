"""Deterministic GFF3 export for parsed entities and localized findings."""

from __future__ import annotations

from urllib.parse import quote

from antismash_review.exporters.track_features import (
    TrackFeature,
    iter_track_features,
    stable_feature_id,
)
from antismash_review.models import Record


def _strand(strand: int | None) -> str:
    return "+" if strand == 1 else "-" if strand == -1 else "."


def _attributes(values: list[tuple[str, str]]) -> str:
    return ";".join(f"{key}={quote(value, safe='-_.:')}" for key, value in values if value)


def _part_id(logical_id: str, part_index: int, part_count: int) -> str:
    return logical_id if part_count == 1 else f"{logical_id}.part{part_index}"


def _rows_for_feature(
    feature: TrackFeature,
    logical_id: str,
) -> list[str]:
    rows: list[str] = []
    parts = feature.location.parts
    for part_index, part in enumerate(parts, start=1):
        # Internal/Biopython coordinates are zero-based half-open [start, end).
        # GFF3 is one-based inclusive, so only start is incremented.
        start = part.start + 1
        end = part.end
        row_id = _part_id(logical_id, part_index, len(parts))
        attributes = [
            ("ID", row_id),
            ("feature_id", logical_id),
            *feature.attributes,
        ]
        if len(parts) > 1:
            attributes.append(("part", f"{part_index}/{len(parts)}"))
        if feature.location.cross_origin:
            attributes.append(("cross_origin", "true"))
        rows.append(
            "\t".join(
                (
                    feature.seqid,
                    "antismash-review",
                    feature.feature_type,
                    str(start),
                    str(end),
                    ".",
                    _strand(part.strand if part.strand is not None else feature.location.strand),
                    ".",
                    _attributes(attributes),
                )
            )
        )
    return rows


def _unique_id(base: str, used: set[str], counts: dict[str, int]) -> str:
    counts[base] = counts.get(base, 0) + 1
    candidate = base if counts[base] == 1 else f"{base}.{counts[base]}"
    while candidate in used:
        counts[base] += 1
        candidate = f"{base}.{counts[base]}"
    used.add(candidate)
    return candidate


def render_gff3(records: list[Record]) -> str:
    used: set[str] = set()
    counts: dict[str, int] = {}
    lines = ["##gff-version 3"]
    for feature in iter_track_features(records):
        base = stable_feature_id(
            feature.record_id,
            feature.feature_type,
            feature.ordinal,
            feature.preferred_id,
        )
        logical_id = _unique_id(base, used, counts)
        lines.extend(_rows_for_feature(feature, logical_id))
    return "\n".join(lines) + "\n"
