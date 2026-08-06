from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from .discovery import discover
from .exporters.json_export import dumps_records
from .exporters.markdown import render_records
from .exporters.tables import render_tsv
from .genbank import GenBankParseError, parse_genbank


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="antismash-review")
    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect = subparsers.add_parser("inspect", help="review antiSMASH GenBank output")
    inspect.add_argument("input", type=Path)
    inspect.add_argument("--format", choices=("markdown", "json", "tsv"), default="markdown")
    inspect.add_argument("--output", type=Path)
    inspect.add_argument("--lenient", action="store_true")
    inspect.add_argument("--recursive", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        manifest = discover(args.input, recursive=args.recursive)
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    paths = manifest.aggregate_genbanks or manifest.region_genbanks
    if not paths:
        if manifest.json_files:
            message = (
                "native antiSMASH JSON input is not supported; "
                "provide an aggregate or region GenBank file"
            )
        else:
            message = "no GenBank input found"
        print(f"error: {message}", file=sys.stderr)
        return 2

    try:
        records = [record for path in paths for record in parse_genbank(path, lenient=args.lenient)]
    except GenBankParseError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.format == "json":
        output = dumps_records(records)
    elif args.format == "tsv":
        output = render_tsv(records)
    else:
        output = render_records(records)
    if args.output:
        output_path = args.output.resolve()
        if output_path in paths:
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
