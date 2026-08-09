from pathlib import Path

from antismash_review.exporters.entity_tables import render_domain_tsv
from antismash_review.exporters.json_export import dumps_records
from antismash_review.exporters.markdown import render_records
from antismash_review.exporters.tables import render_tsv
from antismash_review.genbank import parse_genbank
from antismash_review.review import review_record
from tests.fixtures.build_fixture import write_synthetic_genbank


def test_synthetic_fixture_preserves_hierarchy_raw_features_and_qualifiers(
    tmp_path: Path,
) -> None:
    record = parse_genbank(write_synthetic_genbank(tmp_path / "synthetic.gbk"))[0]

    assert len(record.regions) == 1
    assert len(record.candidate_clusters) == 1
    assert len(record.protoclusters) == 1
    assert len(record.proto_cores) == 1
    assert len(record.genes) == 1  # standalone gene is raw evidence, not a CDS
    assert [domain.tool for domain in record.domains] == [
        "nrps_pks_domains",
        "tigrfam",
    ]
    assert record.domains[0].specificity == [
        "KR activity: inactive",
        "KR stereochemistry: C2",
    ]
    assert record.domains[0].qualifiers["specificity"] == tuple(record.domains[0].specificity)
    assert len(record.modules[0].missing_domain_ids) == 1
    assert len(record.deduplicated_pfam_hits) == 1
    assert {raw.feature_type for raw in record.raw_features} >= {
        "gene",
        "tRNA",
        "rRNA",
        "repeat_region",
    }

    codes = {diagnostic.code for diagnostic in review_record(record)}
    assert "pseudogene_in_cluster" in codes
    assert "unrecognized_feature_type" in codes
    assert "partial_cds_at_edge" in codes
    assert "module_domain_missing" in codes
    assert "orphan_module_locus" in codes

    assert "pseudogene_in_cluster" in dumps_records([record])
    assert "pseudogene_in_cluster" in render_records([record])
    assert "pseudogene_in_cluster" in render_tsv([record])
    assert "KR stereochemistry: C2" in render_domain_tsv([record])
