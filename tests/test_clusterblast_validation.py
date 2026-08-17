from __future__ import annotations

import hashlib
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
from antismash_review.exporters.json_export import dumps_records
from antismash_review.genbank import parse_genbank
from tests.fixtures.build_fixture import write_synthetic_genbank

FIXTURE = Path(__file__).parent / "fixtures" / "clusterblast" / "text" / "contig_1_c1.txt"


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
    mod["general"]["record_id"] = "contig_1"
    mod["general"]["search_type"] = "clusterblast"
    mod["knowncluster"] = {
        "schema_version": 5,
        "record_id": "contig_1",
        "search_type": "knownclusterblast",
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
        "record_id": "contig_1",
        "search_type": "subclusterblast",
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

    # Module schema version != 1 or 2
    doc = _minimal_document()
    _module(doc)["schema_version"] = 3
    with pytest.raises(ClusterBlastParseError, match="Unsupported ClusterBlast module schema"):
        parse_clusterblast_json(_write_document(tmp_path, doc))

    # Record ID mismatch
    doc = _minimal_document()
    _module(doc)["record_id"] = "mismatch"
    with pytest.raises(ClusterBlastParseError, match="does not match"):
        parse_clusterblast_json(_write_document(tmp_path, doc))

    # Result schema version not in {1, 3, 5}
    doc = _minimal_document()
    _module(doc)["general"]["schema_version"] = 6
    with pytest.raises(ClusterBlastParseError, match="Unsupported ClusterBlast general result"):
        parse_clusterblast_json(_write_document(tmp_path, doc))

    # Region number < 1
    doc = _minimal_document()
    _result(doc)["region_number"] = 0
    with pytest.raises(ClusterBlastParseError, match="must be positive"):
        parse_clusterblast_json(_write_document(tmp_path, doc))


def test_json_module_schema_1_is_accepted(tmp_path: Path) -> None:
    doc = _minimal_document()
    _module(doc)["schema_version"] = 1
    results = parse_clusterblast_json(_write_document(tmp_path, doc))
    assert len(results) == 1
    assert results[0].module_schema_version == 1


def test_json_result_schema_1_accepted_similarity_is_none(tmp_path: Path) -> None:
    doc = _minimal_document()
    _module(doc)["schema_version"] = 1
    _module(doc)["general"]["schema_version"] = 1
    del _hit_details(doc)["similarity"]

    results = parse_clusterblast_json(_write_document(tmp_path, doc))
    assert len(results) == 1
    assert results[0].result_schema_version == 1
    assert results[0].rankings[0].similarity is None


def test_json_result_schema_2_accepted_data_version_preserved(tmp_path: Path) -> None:
    doc = _minimal_document()
    _module(doc)["general"]["schema_version"] = 2
    _module(doc)["general"]["data_version"] = "1.0"
    del _hit_details(doc)["similarity"]

    results = parse_clusterblast_json(_write_document(tmp_path, doc))
    assert len(results) == 1
    assert results[0].result_schema_version == 2
    assert results[0].data_version == "1.0"
    assert results[0].rankings[0].similarity is None


def test_json_result_schema_3_accepted_similarity_preserved(tmp_path: Path) -> None:
    doc = _minimal_document()
    _module(doc)["general"]["schema_version"] = 3
    _hit_details(doc)["similarity"] = 42

    results = parse_clusterblast_json(_write_document(tmp_path, doc))
    assert len(results) == 1
    assert results[0].result_schema_version == 3
    assert results[0].rankings[0].similarity == 42


@pytest.mark.parametrize("schema_val", [0, 4, 6])
def test_json_rejects_unsupported_result_schemas(tmp_path: Path, schema_val: int) -> None:
    doc = _minimal_document()
    _module(doc)["general"]["schema_version"] = schema_val
    with pytest.raises(ClusterBlastParseError, match="Unsupported ClusterBlast general result"):
        parse_clusterblast_json(_write_document(tmp_path, doc))


def test_json_data_version_absent_is_none(tmp_path: Path) -> None:
    doc = _minimal_document()
    results = parse_clusterblast_json(_write_document(tmp_path, doc))
    assert results[0].data_version is None


def test_json_data_version_non_string_rejected(tmp_path: Path) -> None:
    doc = _minimal_document()
    _module(doc)["general"]["data_version"] = 42
    with pytest.raises(ClusterBlastParseError, match="data_version"):
        parse_clusterblast_json(_write_document(tmp_path, doc))


def test_json_section_record_id_mismatch_rejected(tmp_path: Path) -> None:
    doc = _minimal_document()
    _module(doc)["general"]["record_id"] = "WRONG"
    with pytest.raises(ClusterBlastParseError, match="record_id"):
        parse_clusterblast_json(_write_document(tmp_path, doc))


def test_json_section_search_type_mismatch_rejected(tmp_path: Path) -> None:
    doc = _minimal_document()
    # The section key is ``general`` but upstream serializes ``clusterblast``.
    _module(doc)["general"]["search_type"] = "general"
    with pytest.raises(ClusterBlastParseError, match="search_type"):
        parse_clusterblast_json(_write_document(tmp_path, doc))


def test_json_plain_query_gene_without_pipes(tmp_path: Path) -> None:
    doc = _minimal_document()
    _hit_details(doc)["pairings"][0][0] = "GENE_PLAIN_123"

    results = parse_clusterblast_json(_write_document(tmp_path, doc))
    assert len(results) == 1
    assert results[0].rankings[0].pairings[0].query_gene == "GENE_PLAIN_123"


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


def test_clusterblast_json_fixtures_across_versions() -> None:
    fixtures_dir = Path(__file__).parent / "fixtures" / "clusterblast"
    min_dir = fixtures_dir / "minimal"

    v6_file = min_dir / "schema1_minimal.json"
    v7_0_file = min_dir / "schema2_minimal.json"
    v7_1_file = min_dir / "schema3_minimal.json"
    v8_file = min_dir / "schema5_minimal.json"
    legacy_file = min_dir / "clusterblast_compat_module_schema1.json"

    v6_res = parse_clusterblast_json(v6_file)[0]
    v7_0_res = parse_clusterblast_json(v7_0_file)[0]
    v7_1_res = parse_clusterblast_json(v7_1_file)[0]
    v8_res = parse_clusterblast_json(v8_file)[0]
    legacy_res = parse_clusterblast_json(legacy_file)[0]

    # Version-specific schemas matching upstream antiSMASH serializers
    assert v6_res.module_schema_version == 2
    assert v6_res.result_schema_version == 1
    assert v6_res.data_version is None
    assert v6_res.rankings[0].similarity is None

    assert v7_0_res.module_schema_version == 2
    assert v7_0_res.result_schema_version == 2
    assert v7_0_res.data_version == "1.0"
    assert v7_0_res.rankings[0].similarity is None

    assert v7_1_res.module_schema_version == 2
    assert v7_1_res.result_schema_version == 3
    assert v7_1_res.data_version == "1.0"
    assert v7_1_res.rankings[0].similarity == 42

    assert v8_res.module_schema_version == 2
    assert v8_res.result_schema_version == 5
    assert v8_res.data_version == "1.0"
    assert v8_res.rankings[0].similarity == 42

    assert legacy_res.module_schema_version == 1
    assert legacy_res.result_schema_version == 1

    # Normalized parity across shared fields
    for res in (v6_res, v7_0_res, v7_1_res, v8_res, legacy_res):
        assert res.record_id == "SYNTH.1"
        assert res.region_number == 1
        assert res.search_type == "clusterblast"
        assert len(res.rankings) == 1
        hit = res.rankings[0]
        assert hit.accession == "BGC0001000"
        assert hit.description == "Synthetic reference cluster"
        assert hit.cluster_type == "NRPS"
        assert hit.num_hits == 1
        assert hit.core_gene_hits == 1
        assert hit.blast_score == 100.0
        assert hit.synteny_score == 2
        assert hit.core_bonus == 1
        assert len(hit.pairings) == 1
        pairing = hit.pairings[0]
        assert pairing.query_gene == "SYN_CDS_1"
        assert pairing.subject_gene == "orf1"
        assert pairing.percent_identity == 80.0
        assert pairing.blast_score == 100.0
        assert pairing.percent_coverage == 90.0
        assert pairing.evalue == 1e-10
        assert pairing.subject_protein_id == "BGC0001000_1"


@pytest.mark.parametrize(
    (
        "filename",
        "result_schema",
        "data_version",
        "similarity",
        "reference_fields",
        "subject_fields",
        "protein_fields",
    ),
    [
        (
            "antismash_6_1_1_clusterblast.json",
            1,
            None,
            None,
            {"accession", "cluster_label", "proteins", "description", "cluster_type", "tags"},
            {
                "name",
                "genecluster",
                "start",
                "end",
                "strand",
                "annotation",
                "perc_ident",
                "blastscore",
                "perc_coverage",
                "evalue",
                "locus_tag",
            },
            {"name", "locus_tag", "location", "strand", "annotations"},
        ),
        (
            "antismash_7_0_1_clusterblast.json",
            2,
            "1.0",
            None,
            {"accession", "cluster_label", "proteins", "description", "cluster_type", "tags"},
            {
                "name",
                "genecluster",
                "start",
                "end",
                "strand",
                "annotation",
                "perc_ident",
                "blastscore",
                "perc_coverage",
                "evalue",
                "locus_tag",
            },
            {"name", "locus_tag", "location", "strand", "annotations"},
        ),
        (
            "antismash_7_1_0_clusterblast.json",
            3,
            "1.0",
            100,
            {"accession", "cluster_label", "proteins", "description", "cluster_type", "tags"},
            {
                "name",
                "genecluster",
                "start",
                "end",
                "strand",
                "annotation",
                "perc_ident",
                "blastscore",
                "perc_coverage",
                "evalue",
                "locus_tag",
            },
            {"name", "locus_tag", "location", "strand", "annotations"},
        ),
        (
            "antismash_8_0_4_clusterblast.json",
            5,
            "1.0",
            100,
            {
                "accession",
                "cluster_label",
                "proteins",
                "description",
                "cluster_type",
                "tags",
                "start",
                "end",
            },
            {
                "name",
                "genecluster",
                "start",
                "end",
                "strand",
                "annotation",
                "perc_ident",
                "blastscore",
                "perc_coverage",
                "evalue",
                "locus_tag",
                "full_name",
            },
            {
                "full_name",
                "name",
                "locus_tag",
                "location",
                "strand",
                "annotations",
                "draw_start",
                "draw_end",
            },
        ),
    ],
)
def test_serializer_reconstructed_golden_fixtures(
    filename: str,
    result_schema: int,
    data_version: str | None,
    similarity: int | None,
    reference_fields: set[str],
    subject_fields: set[str],
    protein_fields: set[str],
) -> None:
    path = Path(__file__).parent / "fixtures" / "clusterblast" / "golden" / filename
    raw = json.loads(path.read_text(encoding="utf-8"))
    section = raw["records"][0]["modules"]["antismash.modules.clusterblast"]["general"]
    region = section["results"][0]
    reference, score = region["ranking"][0]
    subject = score["pairings"][0][2]

    assert section["record_id"] == "SYNTH.1"
    assert section["search_type"] == "clusterblast"
    assert section["schema_version"] == result_schema
    assert section.get("data_version") == data_version
    assert region["prefix"] == "clusterblast"
    assert set(reference) == reference_fields
    assert set(subject) == subject_fields
    assert set(section["proteins"][0]) == protein_fields
    assert score.get("similarity") == similarity

    result = parse_clusterblast_json(path)[0]
    assert result.record_id == "SYNTH.1"
    assert result.module_schema_version == 2
    assert result.result_schema_version == result_schema
    assert result.data_version == data_version
    assert result.rankings[0].similarity == similarity
    assert result.rankings[0].pairings[0].query_gene == "SYN_CDS_1"


def test_golden_fixture_hashes_are_frozen() -> None:
    golden_dir = Path(__file__).parent / "fixtures" / "clusterblast" / "golden"
    expected = {
        "antismash_6_1_1_clusterblast.json": (
            "f750d4f7049db154c147aaf67d69bc3426c70a555ff0ecd2d585b3c16ad63b7d"
        ),
        "antismash_7_0_1_clusterblast.json": (
            "6172e2e82587f3d2a3fd7b5b0291410f819f8e1e0abf79854505782387cd7671"
        ),
        "antismash_7_1_0_clusterblast.json": (
            "27752d74ce4b84bda8e204d12d659455c546b734b76c6990f2ed491257426ecd"
        ),
        "antismash_8_0_4_clusterblast.json": (
            "9971b8bcd4b5b6eb987e2ff5065051eeabc69defe910433368a71d1c9e3aca92"
        ),
    }
    observed = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(golden_dir.glob("*.json"))
    }
    assert observed == expected


def test_data_version_exported_in_record_json_and_clusterblast_tsv(tmp_path: Path) -> None:
    fixture = (
        Path(__file__).parent
        / "fixtures"
        / "clusterblast"
        / "golden"
        / "antismash_7_0_1_clusterblast.json"
    )
    result = parse_clusterblast_json(fixture)[0]
    records = parse_genbank(write_synthetic_genbank(tmp_path / "synthetic.gbk"))
    records[0].clusterblast_results.append(result)

    exported = json.loads(dumps_records(records))
    exported_result = exported["records"][0]["clusterblast_results"][0]
    assert exported_result["module_schema_version"] == 2
    assert exported_result["result_schema_version"] == 2
    assert exported_result["data_version"] == "1.0"

    lines = render_clusterblast_tsv(records).splitlines()
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t"), strict=True))
    assert row["module_schema_version"] == "2"
    assert row["result_schema_version"] == "2"
    assert row["data_version"] == "1.0"
