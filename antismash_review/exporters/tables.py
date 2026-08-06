from __future__ import annotations

import csv
import io

from antismash_review.models import Record
from antismash_review.review import review_record


def render_tsv(records: list[Record]) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    writer.writerow(
        (
            "filename",
            "record_id",
            "region_products",
            "contig_edge",
            "core_genes",
            "total_genes",
            "nrps_pks_domains",
            "all_domains",
            "modules",
            "raw_pfam_hits",
            "deduplicated_pfam_hits",
            "diagnostics",
        )
    )
    for record in records:
        products = [product for region in record.regions for product in region.products]
        edge_values = [region.contig_edge for region in record.regions]
        edge = (
            "true"
            if any(value is True for value in edge_values)
            else "false"
            if edge_values and all(value is False for value in edge_values)
            else ""
        )
        writer.writerow(
            (
                record.source_path.name,
                record.record_id,
                "; ".join(products),
                edge,
                sum(gene.gene_kind == "biosynthetic" for gene in record.genes),
                len(record.genes),
                len(record.nrps_pks_domains),
                len(record.domains),
                len(record.modules),
                len(record.pfam_hits),
                len(record.deduplicated_pfam_hits),
                "; ".join(item.code for item in review_record(record)),
            )
        )
    return output.getvalue()
