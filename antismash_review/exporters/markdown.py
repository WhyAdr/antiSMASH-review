from antismash_review.models import Record
from antismash_review.review import review_record

SEARCH_TYPE_TITLES = {
    "clusterblast": "ClusterBlast",
    "knownclusterblast": "KnownClusterBlast",
    "subclusterblast": "SubClusterBlast",
}


def _text(value: str | None) -> str:
    return value if value else "not reported"


def _escape_cell(value: object | None) -> str:
    if value is None:
        return ""
    text = str(value).replace("\r\n", " ").replace("\n", " ").replace("\t", " ")
    return text.replace("|", r"\|").strip()


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

        if record.clusterblast_results:
            lines.extend(["", "### ClusterBlast", ""])
            sorted_cb = sorted(
                record.clusterblast_results,
                key=lambda r: (r.region_number, r.search_type),
            )
            for res in sorted_cb:
                title = SEARCH_TYPE_TITLES.get(res.search_type, res.search_type)
                parent_dir = res.source_path.parent.name
                if parent_dir in {"clusterblast", "knownclusterblast", "subclusterblast"}:
                    src_str = f"{parent_dir}/{res.source_path.name}"
                else:
                    src_str = res.source_path.name
                lines.extend(
                    [
                        f"#### Region {res.region_number} — {title}",
                        "",
                        f"- Source: {res.source_format} (`{src_str}`)",
                    ]
                )
                if res.total_hits is not None:
                    lines.append(f"- Total database hits: {res.total_hits:,}")
                lines.append(f"- Ranked hits: {len(res.rankings)}")
                lines.append("")

                if res.rankings:
                    header_line = (
                        "| Rank | Accession | Description | Proteins hit | Score | Similarity |"
                    )
                    separator_line = "|---:|---|---|---:|---:|---:|"
                    lines.extend([header_line, separator_line])
                    for hit in res.rankings[:5]:
                        num_hits_str = str(hit.num_hits) if hit.num_hits is not None else ""
                        score_str = f"{hit.blast_score:.1f}" if hit.blast_score is not None else ""
                        sim_str = str(hit.similarity) if hit.similarity is not None else ""
                        lines.append(
                            f"| {hit.rank} | {_escape_cell(hit.accession)} | "
                            f"{_escape_cell(hit.description)} | {num_hits_str} | "
                            f"{score_str} | {sim_str} |"
                        )
                    if len(res.rankings) > 5:
                        lines.extend(["", f"- Showing first 5 of {len(res.rankings)} hits"])
                    lines.append("")

        diagnostics = review_record(record)
        lines.extend(["### Diagnostics", ""])
        if diagnostics:
            lines.extend(
                f"- `{item.severity.value}` `{item.code}`: {item.message}" for item in diagnostics
            )
        else:
            lines.append("- None")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
