from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from .clusterblast import (
    ClusterBlastParseError,
    attach_clusterblast_results,
    merge_clusterblast_results,
    parse_clusterblast_json,
    parse_clusterblast_text,
)
from .compare import compare_records
from .discovery import InputManifest, discover
from .exporters.compare_json import dumps_comparison
from .exporters.compare_markdown import render_comparison
from .exporters.entity_tables import (
    render_clusterblast_tsv,
    render_domain_tsv,
    render_gene_tsv,
)
from .exporters.json_export import dumps_records
from .exporters.markdown import render_records
from .exporters.tables import render_tsv
from .genbank import GenBankParseError, parse_genbank
from .models import Diagnostic, Record, Severity


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="antismash-review")
    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect = subparsers.add_parser("inspect", help="review antiSMASH GenBank output")
    inspect.add_argument("input", type=Path)
    inspect.add_argument(
        "--format",
        choices=("markdown", "json", "tsv", "gene-tsv", "domain-tsv", "clusterblast-tsv"),
        default="markdown",
    )
    inspect.add_argument("--output", type=Path)
    inspect.add_argument("--lenient", action="store_true")
    inspect.add_argument("--recursive", action="store_true")

    compare_cmd = subparsers.add_parser("compare", help="compare two review inputs")
    compare_cmd.add_argument("left", type=Path)
    compare_cmd.add_argument("right", type=Path)
    compare_cmd.add_argument(
        "--match-by",
        choices=("record_id", "record_region", "single_record", "coordinate_overlap"),
        default="record_id",
    )
    compare_cmd.add_argument(
        "--assume-shared-coordinate-system",
        action="store_true",
        help=(
            "confirm that left and right coordinates are directly comparable; "
            "required with --match-by coordinate_overlap"
        ),
    )
    compare_cmd.add_argument(
        "--min-reciprocal-overlap",
        type=float,
        default=0.80,
        metavar="FRACTION",
        help="minimum overlap fraction on both records in coordinate mode (default: 0.80)",
    )
    compare_cmd.add_argument(
        "--format",
        choices=("markdown", "json"),
        default="markdown",
    )
    compare_cmd.add_argument("--output", type=Path)
    compare_cmd.add_argument("--lenient", action="store_true")
    compare_cmd.add_argument("--recursive", action="store_true")

    return parser


def load_review_records(
    manifest: InputManifest,
    *,
    lenient: bool,
) -> tuple[list[Record], set[Path]]:
    paths = manifest.aggregate_genbanks or manifest.region_genbanks
    all_inputs = (
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

    text_results = []
    json_results = []
    try:
        for p in manifest.clusterblast_files:
            text_results.append(parse_clusterblast_text(p, search_type="clusterblast"))
        for p in manifest.knownclusterblast_files:
            text_results.append(parse_clusterblast_text(p, search_type="knownclusterblast"))
        for p in manifest.subclusterblast_files:
            text_results.append(parse_clusterblast_text(p, search_type="subclusterblast"))
        for p in manifest.json_files:
            json_results.extend(parse_clusterblast_json(p))
        if text_results or json_results:
            merged = merge_clusterblast_results(text_results, json_results)
            attach_clusterblast_results(records, merged)
    except ClusterBlastParseError as exc:
        if lenient:
            if records:
                records[0].diagnostics.append(
                    Diagnostic(
                        code="clusterblast_parse_failed",
                        severity=Severity.WARNING,
                        message=str(exc),
                        source=str(manifest.root),
                        record_id=records[0].record_id,
                    )
                )
        else:
            raise

    return records, all_inputs


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "compare":
        if args.match_by == "coordinate_overlap":
            if not args.assume_shared_coordinate_system:
                parser.error(
                    "--match-by coordinate_overlap requires --assume-shared-coordinate-system"
                )
            if not 0 < args.min_reciprocal_overlap <= 1:
                parser.error("--min-reciprocal-overlap must be in the interval (0, 1]")
        elif args.assume_shared_coordinate_system:
            parser.error(
                "--assume-shared-coordinate-system is only valid with --match-by coordinate_overlap"
            )

        try:
            left_manifest = discover(args.left, recursive=args.recursive)
            right_manifest = discover(args.right, recursive=args.recursive)
        except (FileNotFoundError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

        try:
            left_records, left_inputs = load_review_records(left_manifest, lenient=args.lenient)
            right_records, right_inputs = load_review_records(right_manifest, lenient=args.lenient)
        except (GenBankParseError, ClusterBlastParseError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

        try:
            comparison = compare_records(
                left_records,
                right_records,
                left_input=args.left,
                right_input=args.right,
                match_method=args.match_by,
                assume_shared_coordinate_system=args.assume_shared_coordinate_system,
                min_reciprocal_overlap=args.min_reciprocal_overlap,
            )
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

        if args.format == "json":
            output = dumps_comparison(comparison)
        else:
            output = render_comparison(comparison)

        if args.output:
            output_path = args.output.resolve()
            all_inputs = left_inputs | right_inputs
            if output_path in all_inputs:
                print(f"error: refusing to overwrite input file: {output_path}", file=sys.stderr)
                return 2
            try:
                args.output.write_text(output, encoding="utf-8")
            except OSError as exc:
                print(f"error: could not write {args.output}: {exc}", file=sys.stderr)
                return 2
        else:
            print(output, end="")
        return 0

    try:
        manifest = discover(args.input, recursive=args.recursive)
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    try:
        records, all_inputs = load_review_records(manifest, lenient=args.lenient)
    except (GenBankParseError, ClusterBlastParseError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.format == "json":
        output = dumps_records(records)
    elif args.format == "tsv":
        output = render_tsv(records)
    elif args.format == "gene-tsv":
        output = render_gene_tsv(records)
    elif args.format == "domain-tsv":
        output = render_domain_tsv(records)
    elif args.format == "clusterblast-tsv":
        output = render_clusterblast_tsv(records)
    else:
        output = render_records(records)

    if args.output:
        output_path = args.output.resolve()
        if output_path in all_inputs:
            print(f"error: refusing to overwrite input file: {output_path}", file=sys.stderr)
            return 2
        try:
            args.output.write_text(output, encoding="utf-8")
        except OSError as exc:
            print(f"error: could not write {args.output}: {exc}", file=sys.stderr)
            return 2
    else:
        print(output, end="")
    return 0
