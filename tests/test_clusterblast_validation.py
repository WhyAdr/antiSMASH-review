from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from antismash_review.clusterblast import ClusterBlastParseError, parse_clusterblast_json


def _minimal_document() -> dict[str, Any]:
    return {
        "records": [
            {
                "id": "contig_1",
                "modules": {
                    "antismash.modules.clusterblast": {
                        "schema_version": 2,
                        "record_id": "contig_1",
                        "general": {
                            "schema_version": 5,
                            "results": [
                                {
                                    "region_number": 1,
                                    "total_hits": 1,
                                    "ranking": [
                                        [
                                            {
                                                "accession": "ACC1",
                                                "description": "description",
                                                "cluster_type": "NRPS",
                                            },
                                            {
                                                "hits": 1,
                                                "core_gene_hits": 0,
                                                "blast_score": 5.0,
                                                "synteny_score": 1,
                                                "core_bonus": 0,
                                                "similarity": 50,
                                                "pairings": [],
                                            },
                                        ]
                                    ],
                                }
                            ],
                        },
                    }
                },
            }
        ]
    }


def _module(document: dict[str, Any]) -> dict[str, Any]:
    return document["records"][0]["modules"]["antismash.modules.clusterblast"]


def _result(document: dict[str, Any]) -> dict[str, Any]:
    return _module(document)["general"]["results"][0]


def _hit_details(document: dict[str, Any]) -> dict[str, Any]:
    return _result(document)["ranking"][0][1]


def _write_document(tmp_path: Path, document: dict[str, Any]) -> Path:
    path = tmp_path / "clusterblast.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("section", "value"),
    [
        ("module", 2.0),
        ("module", "2"),
        ("module", True),
        ("module", None),
        ("result", 5.0),
        ("result", "5"),
        ("result", True),
        ("result", None),
    ],
)
def test_json_clusterblast_rejects_non_integer_schema_versions(
    tmp_path: Path,
    section: str,
    value: object,
) -> None:
    document = _minimal_document()
    if section == "module":
        _module(document)["schema_version"] = value
    else:
        _module(document)["general"]["schema_version"] = value

    with pytest.raises(ClusterBlastParseError, match=r"schema_version is not an integer"):
        parse_clusterblast_json(_write_document(tmp_path, document))


@pytest.mark.parametrize(
    "field",
    ["total_hits", "hits", "core_gene_hits", "similarity"],
)
def test_json_clusterblast_rejects_negative_counts_and_similarity(
    tmp_path: Path,
    field: str,
) -> None:
    document = _minimal_document()
    if field == "total_hits":
        _result(document)[field] = -1
    else:
        _hit_details(document)[field] = -1

    with pytest.raises(ClusterBlastParseError, match=r"must not be negative"):
        parse_clusterblast_json(_write_document(tmp_path, document))


def test_json_clusterblast_preserves_integer_schema_versions_and_zero_counts(
    tmp_path: Path,
) -> None:
    document = _minimal_document()
    _result(document)["total_hits"] = 0
    _hit_details(document)["hits"] = 0
    _hit_details(document)["core_gene_hits"] = 0
    _hit_details(document)["similarity"] = 0

    parsed = parse_clusterblast_json(_write_document(tmp_path, document))[0]

    assert type(parsed.module_schema_version) is int
    assert type(parsed.result_schema_version) is int
    assert parsed.total_hits == 0
    assert parsed.rankings[0].num_hits == 0
    assert parsed.rankings[0].core_gene_hits == 0
    assert parsed.rankings[0].similarity == 0
