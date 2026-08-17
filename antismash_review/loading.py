"""Canonical discovery, parsing, and ClusterBlast sidecar loading."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .clusterblast import (
    ClusterBlastParseError,
    attach_clusterblast_results,
    merge_clusterblast_results,
    parse_clusterblast_json,
    parse_clusterblast_text,
)
from .discovery import InputManifest
from .genbank import parse_genbank
from .models import ClusterBlastResult, Diagnostic, Record, Severity


@dataclass(slots=True)
class LoadedReviewInput:
    """Records and every file that must be protected from output overwrite."""

    root: Path
    records: list[Record]
    input_paths: set[Path]
    input_diagnostics: list[Diagnostic] = field(default_factory=list)


def load_review_input(
    manifest: InputManifest,
    *,
    lenient: bool = False,
) -> LoadedReviewInput:
    """Load GenBank records and enrich them with supported ClusterBlast sidecars.

    Aggregate GenBank files take precedence over region GenBank files, matching the
    existing discovery contract.  All discovered source and sidecar paths are returned
    so callers can refuse to overwrite an input with an output file.
    """

    paths = manifest.aggregate_genbanks or manifest.region_genbanks
    input_paths = (
        set(manifest.aggregate_genbanks)
        | set(manifest.region_genbanks)
        | set(manifest.json_files)
        | set(manifest.clusterblast_files)
        | set(manifest.knownclusterblast_files)
        | set(manifest.subclusterblast_files)
    )

    if not paths:
        if manifest.json_files:
            message = (
                "native antiSMASH JSON cannot yet provide the review record model; "
                "provide GenBank, optionally in a result directory with JSON used as "
                "ClusterBlast enrichment"
            )
        else:
            message = "no GenBank input found"
        raise ValueError(message)

    records = [record for path in paths for record in parse_genbank(path, lenient=lenient)]

    text_results: list[ClusterBlastResult] = []
    json_results: list[ClusterBlastResult] = []
    input_diagnostics: list[Diagnostic] = []

    def _input_diagnostic(exc: ClusterBlastParseError, source: Path) -> None:
        input_diagnostics.append(
            Diagnostic(
                code="clusterblast_parse_failed",
                severity=Severity.WARNING,
                message=str(exc),
                source=str(source),
                record_id=None,
            )
        )

    for path in manifest.clusterblast_files:
        try:
            text_results.append(parse_clusterblast_text(path, search_type="clusterblast"))
        except ClusterBlastParseError as exc:
            if not lenient:
                raise
            _input_diagnostic(exc, path)

    for path in manifest.knownclusterblast_files:
        try:
            text_results.append(parse_clusterblast_text(path, search_type="knownclusterblast"))
        except ClusterBlastParseError as exc:
            if not lenient:
                raise
            _input_diagnostic(exc, path)

    for path in manifest.subclusterblast_files:
        try:
            text_results.append(parse_clusterblast_text(path, search_type="subclusterblast"))
        except ClusterBlastParseError as exc:
            if not lenient:
                raise
            _input_diagnostic(exc, path)

    for path in manifest.json_files:
        try:
            json_results.extend(parse_clusterblast_json(path))
        except ClusterBlastParseError as exc:
            if not lenient:
                raise
            _input_diagnostic(exc, path)

    if text_results or json_results:
        merge_diagnostics: list[Diagnostic] = []
        attach_diagnostics: list[Diagnostic] = []
        try:
            merged = merge_clusterblast_results(
                text_results, json_results, lenient=lenient, diagnostics=merge_diagnostics
            )
            for diag in merge_diagnostics:
                matches = [r for r in records if diag.record_id in {r.record_id, r.name}]
                if len(matches) == 1:
                    matches[0].diagnostics.append(diag)
                else:
                    input_diagnostics.append(diag)
            attach_clusterblast_results(
                records, merged, lenient=lenient, diagnostics=attach_diagnostics
            )
            input_diagnostics.extend(attach_diagnostics)
        except ClusterBlastParseError as exc:
            if not lenient:
                raise
            _input_diagnostic(exc, manifest.root)

    return LoadedReviewInput(
        root=manifest.root,
        records=records,
        input_paths=input_paths,
        input_diagnostics=input_diagnostics,
    )
