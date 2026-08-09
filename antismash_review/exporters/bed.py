"""Deterministic BED6 export for parsed entities and localized findings."""

from __future__ import annotations

from antismash_review.exporters.gff3 import _unique_id
from antismash_review.exporters.track_features import iter_track_features, stable_feature_id
from antismash_review.models import Record


def _strand(strand: int | None) -> str:
    return "+" if strand == 1 else "-" if strand == -1 else "."


def render_bed(records: list[Record]) -> str:
    used: set[str] = set()
    counts: dict[str, int] = {}
    lines: list[str] = []
    for feature in iter_track_features(records):
        base = stable_feature_id(
            feature.record_id,
            feature.feature_type,
            feature.ordinal,
            feature.preferred_id,
        )
        logical_id = _unique_id(base, used, counts)
        parts = feature.location.parts
        for part_index, part in enumerate(parts, start=1):
            name = logical_id if len(parts) == 1 else f"{logical_id}.part{part_index}"
            lines.append(
                "\t".join(
                    (
                        feature.seqid,
                        str(part.start),
                        str(part.end),
                        f"{feature.feature_type}|{name}",
                        "0",
                        _strand(
                            part.strand if part.strand is not None else feature.location.strand
                        ),
                    )
                )
            )
    return "\n".join(lines) + ("\n" if lines else "")
