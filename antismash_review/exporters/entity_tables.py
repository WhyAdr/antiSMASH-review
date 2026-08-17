from __future__ import annotations

import csv
import io
from collections.abc import Sequence

from antismash_review.models import Record

GENE_COLUMNS = (
    "source_path",
    "source_sha256",
    "filename",
    "record_id",
    "locus_tag",
    "gene",
    "product",
    "gene_kind",
    "start",
    "end",
    "strand",
    "partial",
    "cross_origin",
    "ec_numbers",
    "db_xrefs",
    "region_numbers",
    "candidate_cluster_numbers",
    "protocluster_numbers",
    "proto_core_numbers",
)

DOMAIN_COLUMNS = (
    "source_path",
    "source_sha256",
    "filename",
    "record_id",
    "domain_id",
    "name",
    "tool",
    "is_nrps_pks",
    "locus_tag",
    "start",
    "end",
    "strand",
    "partial",
    "cross_origin",
    "protein_start",
    "protein_end",
    "score",
    "evalue",
    "subtypes",
    "specificity",
)

CLUSTERBLAST_COLUMNS = (
    "record_id",
    "region_number",
    "search_type",
    "result_source",
    "result_source_path",
    "total_hits",
    "ranked_hit_count",
    "rank",
    "accession",
    "description",
    "cluster_type",
    "num_hits",
    "blast_score",
    "similarity",
    "pairing_count",
    "module_schema_version",
    "result_schema_version",
    "data_version",
)


def _list_cell(values: Sequence[object]) -> str:
    return "; ".join(str(value) for value in values)


def _optional(value: object | None) -> object:
    return "" if value is None else value


def _boolean(value: bool) -> str:
    return "true" if value else "false"


def render_gene_tsv(records: list[Record]) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    writer.writerow(GENE_COLUMNS)
    for record in records:
        for gene in record.genes:
            writer.writerow(
                (
                    str(record.source_path),
                    record.source_sha256,
                    record.source_path.name,
                    record.record_id,
                    _optional(gene.locus_tag),
                    _optional(gene.gene),
                    _optional(gene.product),
                    gene.gene_kind,
                    gene.location.start,
                    gene.location.end,
                    _optional(gene.location.strand),
                    _boolean(gene.location.partial),
                    _boolean(gene.location.cross_origin),
                    _list_cell(gene.ec_numbers),
                    _list_cell(gene.db_xrefs),
                    _list_cell(gene.region_numbers),
                    _list_cell(gene.candidate_cluster_numbers),
                    _list_cell(gene.protocluster_numbers),
                    _list_cell(gene.proto_core_numbers),
                )
            )
    return output.getvalue()


def render_domain_tsv(records: list[Record]) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    writer.writerow(DOMAIN_COLUMNS)
    for record in records:
        for domain in record.domains:
            writer.writerow(
                (
                    str(record.source_path),
                    record.source_sha256,
                    record.source_path.name,
                    record.record_id,
                    _optional(domain.domain_id),
                    _optional(domain.name),
                    _optional(domain.tool),
                    _boolean(domain.is_nrps_pks),
                    _optional(domain.locus_tag),
                    domain.location.start,
                    domain.location.end,
                    _optional(domain.location.strand),
                    _boolean(domain.location.partial),
                    _boolean(domain.location.cross_origin),
                    _optional(domain.protein_start),
                    _optional(domain.protein_end),
                    _optional(domain.score),
                    _optional(domain.evalue),
                    _list_cell(domain.subtypes),
                    _list_cell(domain.specificity),
                )
            )
    return output.getvalue()


def render_clusterblast_tsv(records: list[Record]) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    writer.writerow(CLUSTERBLAST_COLUMNS)
    for record in records:
        for result in record.clusterblast_results:
            if not result.rankings:
                writer.writerow(
                    (
                        result.record_id,
                        result.region_number,
                        result.search_type,
                        result.source_format,
                        str(result.source_path),
                        _optional(result.total_hits),
                        0,
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        _optional(result.module_schema_version),
                        _optional(result.result_schema_version),
                        _optional(result.data_version),
                    )
                )
            else:
                for hit in result.rankings:
                    writer.writerow(
                        (
                            result.record_id,
                            result.region_number,
                            result.search_type,
                            result.source_format,
                            str(result.source_path),
                            _optional(result.total_hits),
                            len(result.rankings),
                            hit.rank,
                            hit.accession,
                            hit.description,
                            _optional(hit.cluster_type),
                            _optional(hit.num_hits),
                            _optional(hit.blast_score),
                            _optional(hit.similarity),
                            len(hit.pairings),
                            _optional(result.module_schema_version),
                            _optional(result.result_schema_version),
                            _optional(result.data_version),
                        )
                    )
    return output.getvalue()
