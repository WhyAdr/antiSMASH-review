from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from antismash_review.genbank import GenBankParseError, _antismash_version, parse_genbank
from antismash_review.models import RawFeature, Record
from antismash_review.review import review_record

ROOT = Path(__file__).resolve().parent / "fixtures"


def test_two_records_remain_distinct() -> None:
    records = parse_genbank(ROOT / "multi-record.gb")
    assert [record.record_id for record in records] == ["FIRST.1", "SECOND.1"]
    assert [record.regions[0].products for record in records] == [["first"], ["second"]]


def test_hierarchy_and_raw_features_are_lossless() -> None:
    record = parse_genbank(ROOT / "semantics.gb")[0]
    assert record.regions[0].products == ["hybrid-a", "hybrid-b"]
    assert len(record.candidate_clusters) == 3
    assert len(record.protoclusters) == 2
    assert record.regions[0].rules == ["rule one", "rule two"]
    assert record.candidate_clusters[0].smiles == ["CC"]
    assert any(feature.feature_type == "mystery_feature" for feature in record.raw_features)
    assert record.antismash_version == "8.0.4"


def test_domains_modules_motifs_and_pfam_views() -> None:
    record = parse_genbank(ROOT / "semantics.gb")[0]
    assert len(record.domains) == 2
    assert len(record.nrps_pks_domains) == 1
    assert record.domains[0].subtypes == ["other"]
    assert record.domains[1].subtypes == ["AMP-binding"]

    module = record.modules[0]
    assert module.domain_ids == ["tigr_domain_1", "nrps_domain_2"]
    assert module.complete and module.starter and module.final and module.iterative
    assert module.multi_cds and not module.missing_domain_ids

    assert record.motifs[0].label == "NRPS-A_a10"
    assert record.motifs[0].evalue == pytest.approx(0.63)
    assert record.motifs[1].label is None
    assert record.motifs[1].core_sequence == "corepeptide"
    assert len(record.pfam_hits) == 2
    assert len(record.deduplicated_pfam_hits) == 1


def test_missing_gene_kind_is_unclassified() -> None:
    record = parse_genbank(ROOT / "semantics.gb")[0]
    assert [gene.gene_kind for gene in record.genes] == ["unclassified", "other", "biosynthetic"]


def test_locations_preserve_parts_fuzziness_and_membership() -> None:
    record = parse_genbank(ROOT / "semantics.gb")[0]
    location = record.proto_cores[0].location
    assert location.cross_origin
    assert location.partial
    assert len(location.parts) == 2
    assert record.genes[0].region_numbers == [7]
    assert record.genes[0].candidate_cluster_numbers == [1, 2, 3]


def test_boundary_review_is_evidence_scoped() -> None:
    record = parse_genbank(ROOT / "semantics.gb")[0]
    diagnostics = review_record(record)
    codes = {item.code for item in diagnostics}
    assert "context_reaches_record_edge" in codes
    assert "core_reaches_record_edge" in codes
    assert "orphan_module_locus" in codes  # G4 not in CDS set
    assert not any(item.code == "low_confidence_motif" for item in diagnostics)


def test_strict_parse_does_not_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    from antismash_review import genbank

    def fail(*args: object, **kwargs: object) -> object:
        raise ValueError("synthetic parser failure")

    monkeypatch.setattr(genbank.SeqIO, "parse", fail)
    with pytest.raises(GenBankParseError, match="synthetic parser failure"):
        parse_genbank(ROOT / "multi-record.gb")


def test_lenient_mode_records_adapter_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    from antismash_review import genbank

    original = genbank._adapt_feature

    def fail_on_cds(record: Record, raw: RawFeature) -> None:
        if raw.feature_type == "CDS":
            raise ValueError("synthetic adapter failure")
        original(record, raw)

    monkeypatch.setattr(genbank, "_adapt_feature", fail_on_cds)
    record = parse_genbank(ROOT / "semantics.gb", lenient=True)[0]
    assert len(record.genes) == 0
    assert any(item.code == "feature_adapter_failed" for item in record.diagnostics)


def test_version_comes_only_from_antismash_metadata() -> None:
    source = SeqRecord(Seq("A"), id="version-test")
    source.annotations["structured_comment"] = {
        "Other-Data": {"Version": "wrong"},
        "antiSMASH-Data": {"Version": "8.0.4"},
    }
    assert _antismash_version(source) == "8.0.4"

    source.annotations["structured_comment"] = {"Other-Data": {"Version": "wrong"}}
    source.annotations["comment"] = (
        "##antiSMASH-Data-START##\nVersion :: 7.1.0\n##antiSMASH-Data-END##"
    )
    assert _antismash_version(source) == "7.1.0"


def test_missing_source_uses_public_parse_error(tmp_path: Path) -> None:
    with pytest.raises(GenBankParseError, match="Could not parse"):
        parse_genbank(tmp_path / "missing.gbk")


def test_empty_genbank_raises_error(tmp_path: Path) -> None:
    empty = tmp_path / "empty.gb"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(GenBankParseError, match="No GenBank records found"):
        parse_genbank(empty)


def test_cross_origin_location_parts_and_membership() -> None:
    record = parse_genbank(ROOT / "cross-origin.gb")[0]
    assert record.topology == "circular"
    gene = record.genes[0]
    assert gene.locus_tag == "XORIGIN1"
    assert gene.location.cross_origin
    assert len(gene.location.parts) == 2
    # CDS spans origin: parts should cover 80..100 and 0..20
    part_ranges = [(p.start, p.end) for p in gene.location.parts]
    assert (80, 100) in part_ranges
    assert (0, 20) in part_ranges
    # Gene should still be assigned to its region
    assert gene.region_numbers == [1]


def test_partial_edge_fixture() -> None:
    record = parse_genbank(ROOT / "partial-edge.gb")[0]
    gene = record.genes[0]
    assert gene.locus_tag == "EDGE1"
    assert gene.location.partial
    # Fuzzy start at <1 means start is 0, which is a boundary
    assert gene.location.start == 0
    diagnostics = review_record(record)
    assert any(item.code == "partial_cds_at_edge" for item in diagnostics)
