from __future__ import annotations

import json
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from antismash_review.cli import main
from antismash_review.clustering import domain_jaccard_distances
from antismash_review.cohort import CohortError, build_cohort
from antismash_review.exporters.cohort_json import dumps_cohort
from antismash_review.exporters.cohort_table import (
    render_domain_matrix_tsv,
    render_product_matrix_tsv,
)


def _write_member(
    path: Path,
    *,
    record_id: str,
    products: list[str],
    domains: list[str],
) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    source = SeqRecord(Seq("A" * 400), id=record_id, name=record_id, description="cohort fixture")
    source.annotations["molecule_type"] = "DNA"
    source.annotations["structured_comment"] = {"antiSMASH-Data": {"Version": "8.0.4"}}
    source.features = [
        SeqFeature(
            FeatureLocation(0, 400, strand=1),
            type="source",
            qualifiers={"organism": ["cohort fixture"]},
        ),
        SeqFeature(
            FeatureLocation(0, 400, strand=1),
            type="region",
            qualifiers={"region_number": ["1"], "product": products},
        ),
    ]
    for index, domain_name in enumerate(domains, start=1):
        source.features.append(
            SeqFeature(
                FeatureLocation(index * 10, index * 10 + 8, strand=1),
                type="aSDomain",
                qualifiers={
                    "domain_id": [f"D{index}"],
                    "aSDomain": [domain_name],
                    "aSTool": ["nrps_pks_domains"],
                },
            )
        )
    output = path / "aggregate.gbk"
    SeqIO.write(source, output, "genbank")
    return output


def test_root_mode_is_sorted_and_matrices_have_explicit_units(tmp_path: Path) -> None:
    _write_member(
        tmp_path / "beta",
        record_id="BETA.1",
        products=["T1PKS"],
        domains=["PKS_KS"],
    )
    _write_member(
        tmp_path / "alpha",
        record_id="ALPHA.1",
        products=[" NRPS ", "NRPS"],
        domains=["A", "PCP"],
    )

    result = build_cohort(tmp_path)

    assert [member.name for member in result.members] == ["alpha", "beta"]
    assert result.product_columns == ["nrps", "t1pks"]
    assert result.product_display_labels["nrps"] == "NRPS"
    assert result.product_raw_labels["nrps"] == (" NRPS ", "NRPS")
    assert result.product_matrix == [[1, 0], [0, 1]]
    assert result.domain_columns == ["a", "pcp", "pks_ks"]
    assert result.domain_matrix == [[1, 1, 0], [0, 0, 1]]

    counts = build_cohort(tmp_path, value_mode="count")
    assert counts.product_matrix == [[2, 0], [0, 1]]

    product_tsv = render_product_matrix_tsv(result)
    domain_tsv = render_domain_matrix_tsv(result)
    assert product_tsv.splitlines() == ["sample\tNRPS\tT1PKS", "alpha\t1\t0", "beta\t0\t1"]
    assert domain_tsv.splitlines()[0] == "sample\tA\tPCP\tPKS_KS"


def test_aggregate_genbank_precedence_avoids_duplicate_region_counts(tmp_path: Path) -> None:
    member = tmp_path / "strain"
    _write_member(member, record_id="STRAIN.1", products=["NRPS"], domains=["A"])
    region = member / "contig.region001.gbk"
    region.write_bytes((member / "aggregate.gbk").read_bytes())

    result = build_cohort(tmp_path)

    assert result.product_columns == ["nrps"]
    assert result.product_matrix == [[1]]
    assert region.resolve() in result.input_paths


def test_manifest_preserves_order_and_rejects_duplicate_names(tmp_path: Path) -> None:
    _write_member(tmp_path / "alpha", record_id="ALPHA.1", products=["NRPS"], domains=[])
    _write_member(tmp_path / "beta", record_id="BETA.1", products=["T1PKS"], domains=[])
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("sample\tpath\nbeta\tbeta\nalpha\talpha\n", encoding="utf-8")

    result = build_cohort(manifest=manifest)
    assert [member.name for member in result.members] == ["beta", "alpha"]
    assert manifest.resolve() in result.input_paths

    duplicate = tmp_path / "duplicate.tsv"
    duplicate.write_text("one\talpha\nONE\tbeta\n", encoding="utf-8")
    with pytest.raises(CohortError, match="duplicate cohort sample name"):
        build_cohort(manifest=duplicate)


def test_invalid_member_is_named_by_default_and_reported_when_skipped(tmp_path: Path) -> None:
    _write_member(tmp_path / "valid", record_id="VALID.1", products=["NRPS"], domains=[])
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("valid\tvalid\nbroken\tmissing\n", encoding="utf-8")

    with pytest.raises(CohortError, match="broken.*missing"):
        build_cohort(manifest=manifest)

    result = build_cohort(manifest=manifest, skip_invalid_members=True)
    assert [item.name for item in result.skipped] == ["broken"]
    document = json.loads(dumps_cohort(result))
    assert document["schema_name"] == "antismash-review-cohort"
    assert document["schema_version"] == "0.1.0"
    assert document["skipped"][0]["name"] == "broken"
    assert dumps_cohort(result) == dumps_cohort(result)


def test_cohort_cli_dispatches_matrix_and_refuses_manifest_overwrite(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _write_member(tmp_path / "alpha", record_id="ALPHA.1", products=["NRPS"], domains=[])
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("alpha\talpha\n", encoding="utf-8")

    assert main(["cohort", "--manifest", str(manifest)]) == 0
    assert "sample\tNRPS" in capsys.readouterr().out

    assert (
        main(["cohort", "--manifest", str(manifest), "--format", "json", "--output", str(manifest)])
        == 2
    )
    assert "refusing to overwrite input file" in capsys.readouterr().err


def test_cohort_handles_a_small_47_member_smoke_set(tmp_path: Path) -> None:
    for index in range(47):
        _write_member(
            tmp_path / f"sample_{index:02d}",
            record_id=f"S{index:02d}.1",
            products=["NRPS" if index % 2 else "T1PKS"],
            domains=["A" if index % 2 else "PKS_KS"],
        )

    result = build_cohort(tmp_path)
    assert len(result.members) == 47
    assert len(result.product_matrix) == 47


def test_domain_jaccard_clustering_is_hand_calculable_and_deterministic(tmp_path: Path) -> None:
    _write_member(tmp_path / "gamma", record_id="GAMMA.1", products=["NRPS"], domains=["C"])
    _write_member(tmp_path / "alpha", record_id="ALPHA.1", products=["NRPS"], domains=["A", "B"])
    _write_member(tmp_path / "beta", record_id="BETA.1", products=["NRPS"], domains=["A"])

    result = build_cohort(tmp_path, cluster_by="domain-jaccard")

    assert result.domain_jaccard == [
        [0.0, 0.5, 1.0],
        [0.5, 0.0, 1.0],
        [1.0, 1.0, 0.0],
    ]
    assert result.cluster_order == ["alpha", "beta", "gamma"]
    assert result.cluster_newick == "((alpha:0.25,beta:0.25):0.25,gamma:0.5);"
    assert domain_jaccard_distances(result.members) == result.domain_jaccard

    reversed_manifest = tmp_path / "reversed.tsv"
    reversed_manifest.write_text("gamma\tgamma\nbeta\tbeta\nalpha\talpha\n", encoding="utf-8")
    reversed_result = build_cohort(manifest=reversed_manifest, cluster_by="domain-jaccard")
    assert reversed_result.cluster_order == ["alpha", "beta", "gamma"]
    assert reversed_result.cluster_newick == result.cluster_newick


def test_domain_jaccard_cli_writes_tree_separately(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _write_member(tmp_path / "alpha", record_id="ALPHA.1", products=["NRPS"], domains=["A"])
    _write_member(tmp_path / "beta", record_id="BETA.1", products=["NRPS"], domains=["B"])
    tree = tmp_path / "tree.nwk"

    assert (
        main(
            [
                "cohort",
                str(tmp_path),
                "--format",
                "json",
                "--cluster-by",
                "domain-jaccard",
                "--tree-output",
                str(tree),
            ]
        )
        == 0
    )
    output = json.loads(capsys.readouterr().out)
    assert output["clustering"]["method"] == "average-linkage"
    assert tree.read_text(encoding="utf-8").strip().endswith(";")
