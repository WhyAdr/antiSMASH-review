from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from antismash_review.clusterblast import (
    ClusterBlastParseError,
    attach_clusterblast_results,
    merge_clusterblast_results,
    parse_clusterblast_json,
    parse_clusterblast_text,
)
from antismash_review.exporters.entity_tables import render_clusterblast_tsv
from antismash_review.genbank import parse_genbank
from tests.fixtures.build_fixture import write_synthetic_genbank

FIXTURE = Path(__file__).parent / "fixtures" / "clusterblast" / "contig_1_c1.txt"


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
                                                "pairings": [
                                                    [
                                                        "loc|tag|extra|field|SYN_CDS_1|more",
                                                        1,
                                                        {
                                                            "name": "SUBJ_1",
                                                            "perc_ident": 95.0,
                                                            "blastscore": 100.0,
                                                            "perc_coverage": 90.0,
                                                            "evalue": 1e-20,
                                                            "locus_tag": "SUBJ_LOCUS",
                                                        },
                                                    ]
                                                ],
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


def test_text_clusterblast_parses_pairings_and_provenance() -> None:
    result = parse_clusterblast_text(FIXTURE, search_type="clusterblast")

    assert result.record_id == "SYNTH.1"
    assert result.region_number == 1
    assert result.source_format == "text"
    assert result.rankings[0].pairings[0].query_gene == "SYN_CDS_1"
    assert result.rankings[0].pairings[0].percent_coverage == 80.0


def test_clusterblast_results_merge_text_precedence_and_attach(tmp_path: Path) -> None:
    records = parse_genbank(write_synthetic_genbank(tmp_path / "synthetic.gbk"))
    text_result = parse_clusterblast_text(FIXTURE, search_type="clusterblast")
    merged = merge_clusterblast_results([text_result], [])
    attach_clusterblast_results(records, merged)

    assert len(records[0].clusterblast_results) == 1
    assert records[0].clusterblast_results[0].source_format == "text"
    assert records[0].clusterblast_results[0].source_sha256 == text_result.source_sha256

    tsv_output = render_clusterblast_tsv(records)
    assert "SYNTH.1" in tsv_output
    assert "SYNTH-HIT-1" in tsv_output
    assert records[0].clusterblast_results[0].rankings[0].pairings[0].query_gene == "SYN_CDS_1"


def test_json_clusterblast_all_sections(tmp_path: Path) -> None:
    document = _minimal_document()
    mod = _module(document)
    mod["knowncluster"] = {
        "schema_version": 5,
        "results": [
            {
                "region_number": 1,
                "total_hits": 0,
                "ranking": [],
            }
        ],
    }
    mod["subcluster"] = {
        "schema_version": 5,
        "results": [
            {
                "region_number": 1,
                "total_hits": 0,
                "ranking": [],
            }
        ],
    }
    results = parse_clusterblast_json(_write_document(tmp_path, document))
    assert len(results) == 3
    search_types = {r.search_type for r in results}
    assert search_types == {"clusterblast", "knownclusterblast", "subclusterblast"}


def test_text_clusterblast_parse_errors(tmp_path: Path) -> None:
    # No region in name
    with pytest.raises(ClusterBlastParseError, match="Could not extract region"):
        parse_clusterblast_text(tmp_path / "invalid.txt", search_type="clusterblast")

    # Missing file
    with pytest.raises(ClusterBlastParseError, match="Could not read"):
        parse_clusterblast_text(tmp_path / "contig_c1.txt", search_type="clusterblast")

    # Empty file
    empty_file = tmp_path / "empty_c1.txt"
    empty_file.write_text("", encoding="utf-8")
    with pytest.raises(ClusterBlastParseError, match="Empty ClusterBlast file"):
        parse_clusterblast_text(empty_file, search_type="clusterblast")

    # Invalid header
    bad_header = tmp_path / "bad_c1.txt"
    bad_header.write_text("Invalid header line\n", encoding="utf-8")
    with pytest.raises(ClusterBlastParseError, match="Missing or invalid"):
        parse_clusterblast_text(bad_header, search_type="clusterblast")

    # Missing sections
    no_sections = tmp_path / "nosec_c1.txt"
    no_sections.write_text("ClusterBlast scores for REC1\nSome text\n", encoding="utf-8")
    with pytest.raises(ClusterBlastParseError, match="Significant hits or Details"):
        parse_clusterblast_text(no_sections, search_type="clusterblast")


def test_json_clusterblast_structure_errors(tmp_path: Path) -> None:
    # Not valid JSON
    bad_json = tmp_path / "bad.json"
    bad_json.write_text("{not valid json", encoding="utf-8")
    with pytest.raises(ClusterBlastParseError, match="Could not parse JSON"):
        parse_clusterblast_json(bad_json)

    # Missing records
    with pytest.raises(ClusterBlastParseError, match="Invalid antiSMASH JSON root"):
        parse_clusterblast_json(_write_document(tmp_path, {"no_records": 1}))

    # Module schema version != 2
    doc = _minimal_document()
    _module(doc)["schema_version"] = 3
    with pytest.raises(ClusterBlastParseError, match="Unsupported ClusterBlast module schema"):
        parse_clusterblast_json(_write_document(tmp_path, doc))

    # Record ID mismatch
    doc = _minimal_document()
    _module(doc)["record_id"] = "mismatch"
    with pytest.raises(ClusterBlastParseError, match="does not match"):
        parse_clusterblast_json(_write_document(tmp_path, doc))

    # Result schema version != 5
    doc = _minimal_document()
    _module(doc)["general"]["schema_version"] = 6
    with pytest.raises(ClusterBlastParseError, match="Unsupported ClusterBlast general result"):
        parse_clusterblast_json(_write_document(tmp_path, doc))

    # Region number < 1
    doc = _minimal_document()
    _result(doc)["region_number"] = 0
    with pytest.raises(ClusterBlastParseError, match="must be positive"):
        parse_clusterblast_json(_write_document(tmp_path, doc))


def test_merge_and_attach_errors(tmp_path: Path) -> None:
    text_result = parse_clusterblast_text(FIXTURE, search_type="clusterblast")

    # Duplicate JSON results
    with pytest.raises(ClusterBlastParseError, match="duplicate JSON"):
        merge_clusterblast_results([], [text_result, text_result])

    # Duplicate text results
    with pytest.raises(ClusterBlastParseError, match="duplicate text"):
        merge_clusterblast_results([text_result, text_result], [])

    # Attach error: no candidate records
    records = parse_genbank(write_synthetic_genbank(tmp_path / "synthetic.gbk"))
    mismatch_result = parse_clusterblast_text(FIXTURE, search_type="clusterblast")
    mismatch_result.record_id = "NON_EXISTENT"
    with pytest.raises(ClusterBlastParseError, match="expected one GenBank target"):
        attach_clusterblast_results(records, [mismatch_result])
