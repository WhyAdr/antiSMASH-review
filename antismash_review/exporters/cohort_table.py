"""TSV renderers for cohort product and domain matrices."""

from __future__ import annotations

import csv
import io

from antismash_review.cohort import CohortResult


def _render_matrix(
    result: CohortResult,
    *,
    columns: list[str],
    display_labels: dict[str, str],
    matrix: list[list[int]],
) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample", *(display_labels[column] for column in columns)])
    for member, row in zip(result.members, matrix, strict=True):
        writer.writerow([member.name, *row])
    return output.getvalue()


def render_product_matrix_tsv(result: CohortResult) -> str:
    return _render_matrix(
        result,
        columns=result.product_columns,
        display_labels=result.product_display_labels,
        matrix=result.product_matrix,
    )


def render_domain_matrix_tsv(result: CohortResult) -> str:
    return _render_matrix(
        result,
        columns=result.domain_columns,
        display_labels=result.domain_display_labels,
        matrix=result.domain_matrix,
    )
