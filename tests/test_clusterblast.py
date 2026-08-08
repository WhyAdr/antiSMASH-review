from __future__ import annotations

import csv
import io
import json
import shutil
from pathlib import Path

import pytest

from antismash_review.cli import main
from antismash_review.clusterblast import (
    ClusterBlastParseError,
    attach_clusterblast_results,
    merge_clusterblast_results,
    parse_clusterblast_json,
    parse_clusterblast_text,
)
from antismash_review.exporters.entity_tables import render_clusterblast_tsv
from antismash_review.exporters.markdown import render_records
from antismash_review.models import (
    ClusterBlastHit,
    ClusterBlastResult,
    CollectionFeature,
    Location,
    LocationPart,
    Record,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _minimal_record(record_id: str = "contig_1", **kwargs: object) -> Record:
    defaults: dict[str, object] = {
        "record_id": record_id,
        "name": record_id,
        "description": "test record",
        "length": 10000,
        "molecule_type": "DNA",
        "topology": "linear",
        "source_path": Path("/data/test.gbk"),
        "source_sha256": "abc123sha",
        "antismash_version": "8.0.4",
        "organism": "synthetic",
        "taxonomy": [],
    }
    defaults.update(kwargs)
    return Record(**defaults)  # type: ignore[arg-type]


def _simple_location(start: int, end: int) -> Location:
    part = LocationPart(start=start, end=end, strand=1)
    return Location(
        start=start,
        end=end,
        strand=1,
        parts=(part,),
        cross_origin=False,
        original=f"{start}..{end}",
    )


def test_parse_clusterblast_text_full(tmp_path: Path) -> None:
    txt = tmp_path / "contig_1_c1.txt"
    txt.write_text(
        "ClusterBlast scores for contig_1\n\n"
        "Table of genes:\nZ1 1 10 + g1\n\n"
        "Significant hits:  \n"
        "1. ACC1\tDesc1\n"
        "2. ACC2\tDesc2\n\n"
        "Details:  \n\n"
        ">>\n"
        "1. ACC1\n"
        "Source: Desc1\n"
        "Type: NRPS\n"
        "Number of proteins with BLAST hits to this cluster: 5\n"
        "Cumulative BLAST score: 1234.5\n\n"
        "Table of Blast hits "
        "(query gene, subject gene, %identity, blast score, %coverage, e-value):\n"
        "Q1\tS1\t99.5\t500.0\t100.0\t1.2e-50\n\n"
        ">>\n"
        "2. ACC2\n"
        "Source: Desc2\n"
        "Type: PKS\n"
        "Number of proteins with BLAST hits to this cluster: 3\n"
        "Cumulative BLAST score: 800.0\n\n"
        "Table of Blast hits "
        "(query gene, subject gene, %identity, blast score, %coverage, e-value):\n"
        "Q2\tS2\t80.0\t300.0\t90.0\t2.5e-30\n",
        encoding="utf-8",
    )
    result = parse_clusterblast_text(txt, search_type="clusterblast")
    assert result.record_id == "contig_1"
    assert result.region_number == 1
    assert result.search_type == "clusterblast"
    assert result.total_hits is None
    assert result.source_format == "text"
    assert len(result.rankings) == 2

    hit1 = result.rankings[0]
    assert hit1.rank == 1
    assert hit1.accession == "ACC1"
    assert hit1.description == "Desc1"
    assert hit1.cluster_type == "NRPS"
    assert hit1.num_hits == 5
    assert hit1.blast_score == 1234.5
    assert len(hit1.pairings) == 1
    assert hit1.pairings[0].query_gene == "Q1"
    assert hit1.pairings[0].percent_identity == 99.5


def test_parse_clusterblast_text_empty(tmp_path: Path) -> None:
    txt = tmp_path / "contig_1_c2.txt"
    txt.write_text(
        "ClusterBlast scores for contig_1\n\n"
        "Table of genes:\nZ1 1 10 + g1\n\n"
        "Significant hits:\n\n\n"
        "Details:\n",
        encoding="utf-8",
    )
    result = parse_clusterblast_text(txt, search_type="subclusterblast")
    assert result.region_number == 2
    assert result.search_type == "subclusterblast"
    assert result.rankings == []


def test_parse_clusterblast_text_errors(tmp_path: Path) -> None:
    # Malformed filename without _cN.txt
    bad_name = tmp_path / "invalid.txt"
    bad_name.write_text("ClusterBlast scores for contig_1\nSignificant hits:\nDetails:\n")
    with pytest.raises(ClusterBlastParseError, match="Could not extract region number"):
        parse_clusterblast_text(bad_name, search_type="clusterblast")

    # Non-consecutive ranks
    txt2 = tmp_path / "contig_1_c1.txt"
    txt2.write_text(
        "ClusterBlast scores for contig_1\n\nSignificant hits:\n2. ACC1\tDesc1\n\nDetails:\n",
        encoding="utf-8",
    )
    with pytest.raises(ClusterBlastParseError, match="Non-consecutive rank"):
        parse_clusterblast_text(txt2, search_type="clusterblast")

    # Mismatched details accession
    txt3 = tmp_path / "contig_1_c3.txt"
    txt3.write_text(
        "ClusterBlast scores for contig_1\n\n"
        "Significant hits:\n1. ACC1\tDesc1\n\n"
        "Details:\n\n>>\n1. ACC2\n",
        encoding="utf-8",
    )
    with pytest.raises(ClusterBlastParseError, match="does not match"):
        parse_clusterblast_text(txt3, search_type="clusterblast")

    # Malformed blast hit row
    txt4 = tmp_path / "contig_1_c4.txt"
    txt4.write_text(
        "ClusterBlast scores for contig_1\n\n"
        "Significant hits:\n1. ACC1\tDesc1\n\n"
        "Details:\n\n>>\n1. ACC1\n"
        "Table of Blast hits "
        "(query gene, subject gene, %identity, blast score, %coverage, e-value):\n"
        "Q1\tS1\tinvalid_float\t500.0\t100.0\t1.2e-50\n",
        encoding="utf-8",
    )
    with pytest.raises(ClusterBlastParseError, match="Invalid numeric value"):
        parse_clusterblast_text(txt4, search_type="clusterblast")


def test_parse_clusterblast_json_full(tmp_path: Path) -> None:
    json_file = tmp_path / "results.json"
    doc = {
        "records": [
            {
                "id": "contig_1",
                "modules": {
                    "antismash.modules.clusterblast": {
                        "schema_version": 2,
                        "record_id": "contig_1",
                        "general": {
                            "schema_version": 5,
                            "search_type": "clusterblast",
                            "results": [
                                {
                                    "region_number": 1,
                                    "total_hits": 100,
                                    "ranking": [
                                        [
                                            {
                                                "accession": "NZ_ACC1",
                                                "description": "General hit 1",
                                                "cluster_type": "NRPS",
                                            },
                                            {
                                                "hits": 10,
                                                "core_gene_hits": 2,
                                                "blast_score": 1500.0,
                                                "synteny_score": 12,
                                                "core_bonus": 5,
                                                "similarity": 85,
                                                "pairings": [
                                                    [
                                                        "input|c1|1-100|+|Q_LOCUS|NAD",
                                                        0,
                                                        {
                                                            "name": "S_GENE",
                                                            "perc_ident": 95.0,
                                                            "blastscore": 400.0,
                                                            "perc_coverage": 99.0,
                                                            "evalue": 1e-100,
                                                            "locus_tag": "WP_001",
                                                        },
                                                    ]
                                                ],
                                            },
                                        ]
                                    ],
                                }
                            ],
                        },
                        "knowncluster": {
                            "schema_version": 5,
                            "search_type": "knownclusterblast",
                            "results": [
                                {
                                    "region_number": 1,
                                    "total_hits": 20,
                                    "ranking": [],
                                }
                            ],
                        },
                    }
                },
            }
        ]
    }
    json_file.write_text(json.dumps(doc), encoding="utf-8")
    results = parse_clusterblast_json(json_file)
    assert len(results) == 2
    r1 = [r for r in results if r.search_type == "clusterblast"][0]
    assert r1.record_id == "contig_1"
    assert r1.region_number == 1
    assert r1.total_hits == 100
    assert r1.module_schema_version == 2
    assert r1.result_schema_version == 5
    assert len(r1.rankings) == 1
    hit = r1.rankings[0]
    assert hit.accession == "NZ_ACC1"
    assert hit.num_hits == 10
    assert hit.core_gene_hits == 2
    assert hit.blast_score == 1500.0
    assert hit.synteny_score == 12
    assert hit.core_bonus == 5
    assert hit.similarity == 85
    assert len(hit.pairings) == 1
    p = hit.pairings[0]
    assert p.query_gene == "Q_LOCUS"
    assert p.subject_protein_id == "WP_001"
    assert p.subject_index == 0


def test_parse_clusterblast_json_schema_validation(tmp_path: Path) -> None:
    json_file = tmp_path / "bad_schema.json"
    doc = {
        "records": [
            {
                "id": "contig_1",
                "modules": {
                    "antismash.modules.clusterblast": {
                        "schema_version": 99,
                        "record_id": "contig_1",
                    }
                },
            }
        ]
    }
    json_file.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(
        ClusterBlastParseError, match="Unsupported ClusterBlast module schema version"
    ):
        parse_clusterblast_json(json_file)


def test_merge_clusterblast_precedence_per_key() -> None:
    # Text for (contig_1, 1, clusterblast)
    text_r1 = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="clusterblast",
        total_hits=None,
        rankings=[],
        source_path=Path("cb.txt"),
        source_sha256="sha1",
        source_format="text",
    )
    # JSON for (contig_1, 1, clusterblast) and (contig_1, 1, knownclusterblast)
    json_r1 = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="clusterblast",
        total_hits=50,
        rankings=[],
        source_path=Path("cb.json"),
        source_sha256="sha2",
        source_format="json",
    )
    json_r2 = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="knownclusterblast",
        total_hits=10,
        rankings=[],
        source_path=Path("cb.json"),
        source_sha256="sha2",
        source_format="json",
    )
    merged = merge_clusterblast_results([text_r1], [json_r1, json_r2])
    assert len(merged) == 2
    # The clusterblast key prefers text
    cb_item = [m for m in merged if m.search_type == "clusterblast"][0]
    assert cb_item.source_format == "text"
    # The knownclusterblast key is filled from JSON
    kcb_item = [m for m in merged if m.search_type == "knownclusterblast"][0]
    assert kcb_item.source_format == "json"


def test_attach_clusterblast_results_aggregate_and_region() -> None:
    rec = _minimal_record("contig_1")
    rec.regions.append(
        CollectionFeature(
            feature_type="region",
            number=1,
            location=_simple_location(0, 1000),
            products=["test"],
            references=[],
            kind=None,
            category=None,
            rules=[],
            smiles=[],
            polymer=[],
            core_location=None,
            cutoff=None,
            neighbourhood=None,
            creating_tool=None,
            contig_edge=None,
            qualifiers={},
        )
    )
    res = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="clusterblast",
        total_hits=10,
        rankings=[],
        source_path=Path("cb.txt"),
        source_sha256="sha",
        source_format="text",
    )
    attach_clusterblast_results([rec], [res])
    assert len(rec.clusterblast_results) == 1

    # Unattached result raises error
    res_unattached = ClusterBlastResult(
        record_id="contig_1",
        region_number=99,
        search_type="clusterblast",
        total_hits=10,
        rankings=[],
        source_path=Path("cb.txt"),
        source_sha256="sha",
        source_format="text",
    )
    with pytest.raises(ClusterBlastParseError, match="expected one GenBank target"):
        attach_clusterblast_results([rec], [res_unattached])


def test_clusterblast_tsv_and_markdown_empty_results() -> None:
    rec = _minimal_record("contig_1")
    res = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="subclusterblast",
        total_hits=None,
        rankings=[],
        source_path=Path("subclusterblast/contig_1_c1.txt"),
        source_sha256="sha",
        source_format="text",
    )
    rec.clusterblast_results.append(res)
    tsv = render_clusterblast_tsv([rec])
    reader = csv.reader(io.StringIO(tsv), delimiter="\t")
    rows = list(reader)
    assert len(rows) == 2
    assert rows[1][0] == "contig_1"
    assert rows[1][1] == "1"
    assert rows[1][2] == "subclusterblast"
    assert rows[1][6] == "0"  # ranked_hit_count
    assert rows[1][7] == ""  # rank

    md = render_records([rec])
    assert "### ClusterBlast" in md
    assert "#### Region 1 — SubClusterBlast" in md
    assert "- Ranked hits: 0" in md


def test_clusterblast_markdown_with_hits_and_truncation() -> None:
    rec = _minimal_record("contig_1")
    hits = [
        ClusterBlastHit(
            rank=i,
            accession=f"ACC_{i}",
            description=f"Desc\t{i}|pipe\nnewline",
            cluster_type="NRPS",
            num_hits=5,
            core_gene_hits=1,
            blast_score=100.0 * i,
            synteny_score=i,
            core_bonus=0,
            similarity=90,
            pairings=[],
        )
        for i in range(1, 8)
    ]
    res = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="clusterblast",
        total_hits=1000,
        rankings=hits,
        source_path=Path("clusterblast/contig_1_c1.txt"),
        source_sha256="sha",
        source_format="text",
    )
    rec.clusterblast_results.append(res)
    md = render_records([rec])
    assert "- Total database hits: 1,000" in md
    assert "- Ranked hits: 7" in md
    assert "- Showing first 5 of 7 hits" in md
    assert "ACC_1" in md
    assert "ACC_5" in md
    assert "ACC_6" not in md


def test_merge_clusterblast_duplicate_errors() -> None:
    r1 = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="clusterblast",
        total_hits=None,
        rankings=[],
        source_path=Path("cb1.txt"),
        source_sha256="sha",
        source_format="text",
    )
    r2 = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="clusterblast",
        total_hits=None,
        rankings=[],
        source_path=Path("cb2.txt"),
        source_sha256="sha",
        source_format="text",
    )
    with pytest.raises(ClusterBlastParseError, match="duplicate text ClusterBlast result"):
        merge_clusterblast_results([r1, r2], [])

    jr1 = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="clusterblast",
        total_hits=10,
        rankings=[],
        source_path=Path("cb1.json"),
        source_sha256="sha",
        source_format="json",
    )
    jr2 = ClusterBlastResult(
        record_id="contig_1",
        region_number=1,
        search_type="clusterblast",
        total_hits=10,
        rankings=[],
        source_path=Path("cb2.json"),
        source_sha256="sha",
        source_format="json",
    )
    with pytest.raises(ClusterBlastParseError, match="duplicate JSON ClusterBlast result"):
        merge_clusterblast_results([], [jr1, jr2])


def test_cli_lenient_clusterblast_diagnostic(tmp_path: Path, capsys: object) -> None:
    # Directory with semantics.gb and a corrupt clusterblast text file
    shutil.copy(FIXTURES / "semantics.gb", tmp_path / "semantics.gb")
    cb_dir = tmp_path / "clusterblast"
    cb_dir.mkdir()
    (cb_dir / "contig_1_c1.txt").write_text("CORRUPT TEXT FILE", encoding="utf-8")

    # In strict mode, fails with status 2
    assert main(["inspect", str(tmp_path)]) == 2
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert "error:" in captured.err

    # In lenient mode, succeeds and records diagnostic
    assert main(["inspect", str(tmp_path), "--lenient", "--format", "json"]) == 0
    captured_lenient = capsys.readouterr()  # type: ignore[union-attr]
    doc = json.loads(captured_lenient.out)
    diag_codes = [d["code"] for d in doc["diagnostics"]]
    assert "clusterblast_parse_failed" in diag_codes
