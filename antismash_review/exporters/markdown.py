from __future__ import annotations

from antismash_review.models import Record
from antismash_review.review import review_record


def _text(value: str | None) -> str:
    return value if value else "not reported"


def render_records(records: list[Record]) -> str:
    lines: list[str] = ["# antiSMASH review", ""]
    for record in records:
        lines.extend(
            [
                f"## `{record.record_id}`",
                "",
                f"- Source: `{record.source_path}`",
                f"- antiSMASH: {_text(record.antismash_version)}",
                f"- Length: {record.length:,} bp",
                f"- Regions: {len(record.regions)}",
                f"- CDS features: {len(record.genes)}",
                f"- NRPS/PKS domains: {len(record.nrps_pks_domains)} "
                f"({len(record.domains)} total antiSMASH domains)",
                f"- Modules: {len(record.modules)}",
                f"- Pfam hits: {len(record.pfam_hits)} raw; "
                f"{len(record.deduplicated_pfam_hits)} deduplicated",
                "",
                "### Products",
                "",
            ]
        )
        products = [product for region in record.regions for product in region.products]
        lines.extend(f"- {product}" for product in products)
        if not products:
            lines.append("- No region product reported")

        diagnostics = review_record(record)
        lines.extend(["", "### Diagnostics", ""])
        if diagnostics:
            lines.extend(
                f"- `{item.severity.value}` `{item.code}`: {item.message}" for item in diagnostics
            )
        else:
            lines.append("- None")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
