from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from .clusterblast import ClusterBlastParseError
from .cohort import CohortError, build_cohort
from .compare import compare_records
from .discovery import InputManifest, discover
from .exporters.architecture_json import dumps_architecture
from .exporters.assemblyline_json import dumps_assembly_lines
from .exporters.assemblyline_markdown import render_assemblyline_markdown
from .exporters.assemblyline_table import render_assemblyline_tsv
from .exporters.bed import render_bed
from .exporters.cohort_json import dumps_cohort
from .exporters.cohort_table import render_domain_matrix_tsv, render_product_matrix_tsv
from .exporters.compare_json import dumps_comparison
from .exporters.compare_markdown import render_comparison
from .exporters.entity_tables import (
    render_clusterblast_tsv,
    render_domain_tsv,
    render_gene_tsv,
)
from .exporters.gff3 import render_gff3
from .exporters.json_export import dumps_records
from .exporters.markdown import render_records
from .exporters.provenance import dumps_provenance, render_provenance_tsv
from .exporters.tables import render_tsv
from .genbank import GenBankParseError
from .loading import load_review_input
from .models import Record


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="antismash-review")
    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect = subparsers.add_parser("inspect", help="review antiSMASH GenBank output")
    inspect.add_argument("input", type=Path)
    inspect.add_argument(
        "--format",
        choices=(
            "markdown",
            "json",
            "tsv",
            "gene-tsv",
            "domain-tsv",
            "clusterblast-tsv",
            "assemblyline-tsv",
            "assemblyline-json",
            "assemblyline-markdown",
            "architecture-json",
            "gff3",
            "bed",
            "provenance-json",
            "provenance-tsv",
        ),
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

    cohort_cmd = subparsers.add_parser("cohort", help="build deterministic cohort matrices")
    cohort_cmd.add_argument("root", nargs="?", type=Path)
    cohort_cmd.add_argument("--manifest", type=Path)
    cohort_cmd.add_argument(
        "--format",
        choices=("product-matrix-tsv", "domain-matrix-tsv", "json"),
        default="product-matrix-tsv",
    )
    cohort_cmd.add_argument("--value", choices=("binary", "count"), default="binary")
    cohort_cmd.add_argument("--cluster-by", choices=("none", "domain-jaccard"), default="none")
    cohort_cmd.add_argument("--tree-output", type=Path)
    cohort_cmd.add_argument("--output", type=Path)
    cohort_cmd.add_argument("--lenient", action="store_true")
    cohort_cmd.add_argument("--skip-invalid-members", action="store_true")

    return parser


def load_review_records(
    manifest: InputManifest,
    *,
    lenient: bool = False,
) -> tuple[list[Record], set[Path]]:
    """Compatibility wrapper for the pre-Phase-0 loader API."""

    loaded = load_review_input(manifest, lenient=lenient)
    return loaded.records, loaded.input_paths


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "cohort":
        try:
            cohort = build_cohort(
                args.root,
                manifest=args.manifest,
                value_mode=args.value,
                lenient=args.lenient,
                skip_invalid_members=args.skip_invalid_members,
                cluster_by=args.cluster_by,
            )
        except CohortError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        if args.tree_output and cohort.cluster_newick is None:
            print("error: --tree-output requires --cluster-by domain-jaccard", file=sys.stderr)
            return 2

        if args.format == "product-matrix-tsv":
            output = render_product_matrix_tsv(cohort)
        elif args.format == "domain-matrix-tsv":
            output = render_domain_matrix_tsv(cohort)
        else:
            output = dumps_cohort(cohort)

        if args.output:
            output_path = args.output.resolve()
            if output_path in cohort.input_paths:
                print(f"error: refusing to overwrite input file: {output_path}", file=sys.stderr)
                return 2
            try:
                args.output.write_text(output, encoding="utf-8")
            except OSError as exc:
                print(f"error: could not write {args.output}: {exc}", file=sys.stderr)
                return 2
        else:
            print(output, end="")
        if args.tree_output:
            tree_output = cohort.cluster_newick
            assert tree_output is not None
            tree_path = args.tree_output.resolve()
            if tree_path in cohort.input_paths:
                print(f"error: refusing to overwrite input file: {tree_path}", file=sys.stderr)
                return 2
            try:
                args.tree_output.write_text(tree_output + "\n", encoding="utf-8")
            except OSError as exc:
                print(f"error: could not write {args.tree_output}: {exc}", file=sys.stderr)
                return 2
        return 0

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
            left_loaded = load_review_input(left_manifest, lenient=args.lenient)
            right_loaded = load_review_input(right_manifest, lenient=args.lenient)
        except (GenBankParseError, ClusterBlastParseError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

        try:
            comparison = compare_records(
                left_loaded.records,
                right_loaded.records,
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
            all_inputs = left_loaded.input_paths | right_loaded.input_paths
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
        loaded = load_review_input(manifest, lenient=args.lenient)
    except (GenBankParseError, ClusterBlastParseError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.format == "json":
        output = dumps_records(loaded.records, input_diagnostics=loaded.input_diagnostics)
    elif args.format == "tsv":
        output = render_tsv(loaded.records)
    elif args.format == "gene-tsv":
        output = render_gene_tsv(loaded.records)
    elif args.format == "domain-tsv":
        output = render_domain_tsv(loaded.records)
    elif args.format == "clusterblast-tsv":
        output = render_clusterblast_tsv(loaded.records)
    elif args.format == "assemblyline-tsv":
        output = render_assemblyline_tsv(loaded.records)
    elif args.format == "assemblyline-json":
        output = dumps_assembly_lines(loaded.records)
    elif args.format == "assemblyline-markdown":
        output = render_assemblyline_markdown(loaded.records)
    elif args.format == "architecture-json":
        output = dumps_architecture(loaded.records)
    elif args.format == "gff3":
        output = render_gff3(loaded.records)
    elif args.format == "bed":
        output = render_bed(loaded.records)
    elif args.format == "provenance-json":
        output = dumps_provenance(loaded.records)
    elif args.format == "provenance-tsv":
        output = render_provenance_tsv(loaded.records)
    else:
        output = render_records(loaded.records, input_diagnostics=loaded.input_diagnostics)

    if args.output:
        output_path = args.output.resolve()
        if output_path in loaded.input_paths:
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
