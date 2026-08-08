from __future__ import annotations

from antismash_review.compare import ComparisonResult, IntergenicSummary


def _delta_str(left: int, right: int) -> str:
    diff = right - left
    if diff > 0:
        return f"+{diff}"
    if diff < 0:
        return f"{diff}"
    return "0"


def _fmt_intergenic(summary: IntergenicSummary) -> str:
    if summary.gap_count == 0:
        return "no intergenic gaps"
    mean_s = f"{summary.mean_bp:.1f} bp" if summary.mean_bp is not None else "n/a"
    med_s = f"{summary.median_bp:.1f} bp" if summary.median_bp is not None else "n/a"
    max_s = f"{summary.max_bp:,} bp" if summary.max_bp is not None else "n/a"
    wrap_s = " (circular wrap included)" if summary.circular_wrap_included else ""
    return (
        f"{summary.gap_count} gaps, {summary.total_bp:,} bp total, "
        f"mean {mean_s}, median {med_s}, max {max_s}{wrap_s}"
    )


def render_comparison(result: ComparisonResult) -> str:
    lines: list[str] = [
        "# antiSMASH comparative review",
        "",
        f"- Left input: `{result.left_input}`",
        f"- Right input: `{result.right_input}`",
        f"- Match method: `{result.match_method}`",
    ]

    if result.match_method == "coordinate_overlap":
        lines.append("- Shared coordinate system: asserted")
        if result.min_reciprocal_overlap is not None:
            lines.append(f"- Reciprocal overlap threshold: {result.min_reciprocal_overlap:.0%}")
    elif result.match_method == "single_record":
        lines.append("- Note: single_record is an explicit user-requested pairing")

    lines.append(f"- Matched pairs: {len(result.matched)}")
    lines.append("")

    for comp in result.matched:
        lines.extend(
            [
                f"## Pair: `{comp.match_key}`",
                "",
                f"- Left ID: `{comp.left_record_id}`",
                f"- Right ID: `{comp.right_record_id}`",
            ]
        )

        if comp.coordinate_evidence:
            ev = comp.coordinate_evidence
            lines.extend(
                [
                    f"- Overlap span: {ev.overlap_bp:,} bp",
                    f"- Left coverage: {ev.left_overlap_fraction:.1%} of {ev.left_span_bp:,} bp",
                    f"- Right coverage: {ev.right_overlap_fraction:.1%} of {ev.right_span_bp:,} bp",
                ]
            )

        r_delta = _delta_str(comp.left_region_count, comp.right_region_count)
        g_delta = _delta_str(comp.left_gene_count, comp.right_gene_count)
        d_delta = _delta_str(comp.left_domain_count, comp.right_domain_count)
        np_delta = _delta_str(comp.left_nrps_pks_count, comp.right_nrps_pks_count)
        lines.extend(
            [
                "",
                "### Feature counts",
                "",
                "| Metric | Left | Right | Delta |",
                "|---|---:|---:|---:|",
                f"| Regions | {comp.left_region_count} | {comp.right_region_count} | {r_delta} |",
                f"| Genes | {comp.left_gene_count} | {comp.right_gene_count} | {g_delta} |",
                f"| Domains | {comp.left_domain_count} | {comp.right_domain_count} | {d_delta} |",
                (
                    f"| NRPS/PKS domains | {comp.left_nrps_pks_count} | "
                    f"{comp.right_nrps_pks_count} | {np_delta} |"
                ),
                "",
                "### Products",
                "",
            ]
        )

        if comp.gained_products:
            lines.append(f"- Gained products: {', '.join(comp.gained_products)}")
        else:
            lines.append("- Gained products: None")

        if comp.lost_products:
            lines.append(f"- Lost products: {', '.join(comp.lost_products)}")
        else:
            lines.append("- Lost products: None")

        lines.extend(["", "### Diagnostics", ""])
        if comp.new_diagnostics:
            lines.append("- New diagnostics:")
            for d in comp.new_diagnostics:
                lines.append(f"  - `{d.severity}` `{d.code}`: {d.message}")
        else:
            lines.append("- New diagnostics: None")

        if comp.resolved_diagnostics:
            lines.append("- Resolved diagnostics:")
            for d in comp.resolved_diagnostics:
                lines.append(f"  - `{d.severity}` `{d.code}`: {d.message}")
        else:
            lines.append("- Resolved diagnostics: None")

        lines.extend(
            [
                "",
                "### Intergenic structure",
                "",
                f"- Left: {_fmt_intergenic(comp.left_intergenic)}",
                f"- Right: {_fmt_intergenic(comp.right_intergenic)}",
                "",
            ]
        )

    if result.unmatched_left:
        lines.extend(["## Unmatched left records", ""])
        for uid in result.unmatched_left:
            lines.append(f"- `{uid}`")
        lines.append("")

    if result.unmatched_right:
        lines.extend(["## Unmatched right records", ""])
        for uid in result.unmatched_right:
            lines.append(f"- `{uid}`")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"
