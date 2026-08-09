from pathlib import Path

from antismash_review.discovery import discover


def test_discover_classifies_genbanks_and_prefers_aggregate(tmp_path: Path) -> None:
    (tmp_path / "aggregate.gbk").write_text("", encoding="utf-8")
    (tmp_path / "contig.region001.gbk").write_text("", encoding="utf-8")
    manifest = discover(tmp_path)

    assert manifest.aggregate_genbanks == (tmp_path / "aggregate.gbk",)
    assert manifest.region_genbanks == (tmp_path / "contig.region001.gbk",)


def test_discover_sidecars_only_reads_canonical_top_level_text(tmp_path: Path) -> None:
    for name in ("clusterblast", "knownclusterblast", "subclusterblast"):
        sidecar = tmp_path / name
        sidecar.mkdir()
        direct = sidecar / "contig_1_c1.txt"
        direct.write_text("direct", encoding="utf-8")
        nested = sidecar / "region1"
        nested.mkdir()
        (nested / "details.txt").write_text("html-sidecar placeholder", encoding="utf-8")

    manifest = discover(tmp_path)

    assert manifest.clusterblast_files == (tmp_path / "clusterblast" / "contig_1_c1.txt",)
    assert manifest.knownclusterblast_files == (tmp_path / "knownclusterblast" / "contig_1_c1.txt",)
    assert manifest.subclusterblast_files == (tmp_path / "subclusterblast" / "contig_1_c1.txt",)
    assert all("details.txt" not in str(path) for path in manifest.knownclusterblast_files)
