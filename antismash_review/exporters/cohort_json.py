"""Versioned JSON export for cohort matrices."""

from __future__ import annotations

import json

from antismash_review._version import __version__
from antismash_review.cohort import CohortResult
from antismash_review.exporters.provenance import dumps_provenance
from antismash_review.schema import COHORT_SCHEMA_NAME, COHORT_SCHEMA_VERSION


def _columns(
    keys: list[str],
    display: dict[str, str],
    raw_labels: dict[str, tuple[str, ...]],
) -> list[dict[str, object]]:
    return [
        {"key": key, "display": display[key], "raw_labels": list(raw_labels[key])} for key in keys
    ]


def dumps_cohort(result: CohortResult) -> str:
    members = []
    for member in result.members:
        provenance = json.loads(dumps_provenance(member.records))
        members.append(
            {
                "name": member.name,
                "input_path": str(member.input_path),
                "records": [record.record_id for record in member.records],
                "inputs": provenance["inputs"],
                "product_counts": dict(sorted(member.product_counts.items())),
                "domain_counts": dict(sorted(member.domain_counts.items())),
            }
        )

    document = {
        "schema_name": COHORT_SCHEMA_NAME,
        "schema_version": COHORT_SCHEMA_VERSION,
        "review_tool": {"name": "antismash-review", "version": __version__},
        "root": str(result.root),
        "manifest": str(result.manifest_path) if result.manifest_path is not None else None,
        "value_mode": result.value_mode,
        "normalization": {
            "labels": "Unicode NFKC, strip, casefold",
            "products": "region product labels",
            "domains": "aSDomain name, falling back to domain_id when unnamed",
        },
        "aggregation": {
            "product_counts": "regions across all records in each member",
            "domain_counts": "adapted aSDomain features across all records in each member",
            "binary_values": "member-level presence",
        },
        "columns": {
            "products": _columns(
                result.product_columns,
                result.product_display_labels,
                result.product_raw_labels,
            ),
            "domains": _columns(
                result.domain_columns,
                result.domain_display_labels,
                result.domain_raw_labels,
            ),
        },
        "members": members,
        "matrices": {
            "products": {
                "rows": [member.name for member in result.members],
                "values": result.product_matrix,
            },
            "domains": {
                "rows": [member.name for member in result.members],
                "values": result.domain_matrix,
            },
        },
        "similarity": {
            "metric": "jaccard",
            "domain_values": "binary domain presence",
            "matrix": result.domain_jaccard,
        },
        "clustering": {
            "method": "average-linkage" if result.cluster_order is not None else None,
            "leaf_order": result.cluster_order,
            "newick": result.cluster_newick,
        },
        "skipped": [
            {
                "name": skipped.name,
                "input_path": str(skipped.input_path),
                "error": skipped.error,
            }
            for skipped in result.skipped
        ],
    }
    return json.dumps(document, indent=2, sort_keys=True) + "\n"
