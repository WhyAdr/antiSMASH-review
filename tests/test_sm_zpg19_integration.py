from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from antismash_review.cli import main
from antismash_review.clusterblast import (
    attach_clusterblast_results,
    merge_clusterblast_results,
    parse_clusterblast_json,
    parse_clusterblast_text,
)
from antismash_review.discovery import discover
from antismash_review.genbank import parse_genbank


def test_sm_zpg19_genbank_baseline(sm_zpg19_dir: Path) -> None:
    manifest = discover(sm_zpg19_dir)
    assert len(manifest.aggregate_genbanks) == 1
    gbk_path = manifest.aggregate_genbanks[0]
    records = parse_genbank(gbk_path)
    assert len(records) == 1
    record = records[0]

    assert record.record_id == "contig_1"
    assert record.length == 5269270
    assert record.antismash_version == "8.0.4"
    assert len(record.regions) == 10
    assert len(record.genes) == 4877
    assert len(record.domains) == 178
    assert len(record.nrps_pks_domains) == 56
    assert len(record.modules) == 10
    assert len(record.motifs) == 109
    assert len(record.pfam_hits) == 6967
    assert len(record.deduplicated_pfam_hits) == 6548


def test_sm_zpg19_clusterblast_text_baseline(sm_zpg19_dir: Path) -> None:
    manifest = discover(sm_zpg19_dir)
    assert len(manifest.clusterblast_files) == 10
    assert len(manifest.knownclusterblast_files) == 10
    assert len(manifest.subclusterblast_files) == 10

    cb_results = [
        parse_clusterblast_text(p, search_type="clusterblast") for p in manifest.clusterblast_files
    ]
    kcb_results = [
        parse_clusterblast_text(p, search_type="knownclusterblast")
        for p in manifest.knownclusterblast_files
    ]
    scb_results = [
        parse_clusterblast_text(p, search_type="subclusterblast")
        for p in manifest.subclusterblast_files
    ]

    assert sum(len(r.rankings) for r in cb_results) == 500
    assert sum(len(r.rankings) for r in kcb_results) == 63
    assert sum(len(r.rankings) for r in scb_results) == 0


def test_sm_zpg19_clusterblast_json_baseline(sm_zpg19_dir: Path) -> None:
    json_path = sm_zpg19_dir / "SM-ZPG19_antismash.json"
    if not json_path.exists():
        pytest.skip("SM-ZPG19 JSON file missing")
    json_results = parse_clusterblast_json(json_path)
    assert len(json_results) == 30

    cb = [r for r in json_results if r.search_type == "clusterblast"]
    kcb = [r for r in json_results if r.search_type == "knownclusterblast"]
    scb = [r for r in json_results if r.search_type == "subclusterblast"]

    assert len(cb) == 10
    assert len(kcb) == 10
    assert len(scb) == 10

    assert sum(r.total_hits or 0 for r in cb) == 48710
    assert sum(r.total_hits or 0 for r in kcb) == 63
    assert sum(r.total_hits or 0 for r in scb) == 0

    assert sum(len(r.rankings) for r in cb) == 500
    assert sum(len(r.rankings) for r in kcb) == 63
    assert sum(len(r.rankings) for r in scb) == 0


def test_sm_zpg19_text_json_parity(sm_zpg19_dir: Path) -> None:
    manifest = discover(sm_zpg19_dir)
    json_path = sm_zpg19_dir / "SM-ZPG19_antismash.json"
    if not json_path.exists():
        pytest.skip("SM-ZPG19 JSON file missing")

    text_results = (
        [
            parse_clusterblast_text(p, search_type="clusterblast")
            for p in manifest.clusterblast_files
        ]
        + [
            parse_clusterblast_text(p, search_type="knownclusterblast")
            for p in manifest.knownclusterblast_files
        ]
        + [
            parse_clusterblast_text(p, search_type="subclusterblast")
            for p in manifest.subclusterblast_files
        ]
    )
    json_results = parse_clusterblast_json(json_path)

    text_by_key = {(r.region_number, r.search_type): r for r in text_results}
    json_by_key = {(r.region_number, r.search_type): r for r in json_results}

    assert len(text_by_key) == 30
    assert len(json_by_key) == 30

    for key, text_res in text_by_key.items():
        json_res = json_by_key[key]
        assert len(text_res.rankings) == len(json_res.rankings)
        for t_hit, j_hit in zip(text_res.rankings, json_res.rankings, strict=True):
            assert t_hit.rank == j_hit.rank
            assert t_hit.accession == j_hit.accession
            assert t_hit.description == j_hit.description


def test_sm_zpg19_json_fallback_when_no_text(sm_zpg19_dir: Path, tmp_path: Path) -> None:
    # Copy aggregate GBK and JSON to tmp_path without text folders
    manifest = discover(sm_zpg19_dir)
    gbk_src = manifest.aggregate_genbanks[0]
    json_src = sm_zpg19_dir / "SM-ZPG19_antismash.json"
    if not json_src.exists():
        pytest.skip("SM-ZPG19 JSON file missing")

    shutil.copy(gbk_src, tmp_path / gbk_src.name)
    shutil.copy(json_src, tmp_path / json_src.name)

    temp_manifest = discover(tmp_path)
    assert len(temp_manifest.clusterblast_files) == 0
    assert len(temp_manifest.json_files) == 1

    records = parse_genbank(temp_manifest.aggregate_genbanks[0])
    json_results = parse_clusterblast_json(temp_manifest.json_files[0])
    merged = merge_clusterblast_results([], json_results)
    attach_clusterblast_results(records, merged)

    assert len(records[0].clusterblast_results) == 30
    assert all(r.source_format == "json" for r in records[0].clusterblast_results)


def test_sm_zpg19_json_alone_cli_status_2(
    sm_zpg19_dir: Path, tmp_path: Path, capsys: object
) -> None:
    json_src = sm_zpg19_dir / "SM-ZPG19_antismash.json"
    if not json_src.exists():
        pytest.skip("SM-ZPG19 JSON file missing")

    shutil.copy(json_src, tmp_path / json_src.name)
    result = main(["inspect", str(tmp_path)])
    captured = capsys.readouterr()  # type: ignore[union-attr]
    assert result == 2
    assert "native antiSMASH JSON" in captured.err
