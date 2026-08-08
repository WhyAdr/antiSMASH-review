from pathlib import Path

from antismash_review.genbank import parse_genbank
from antismash_review.models import (
    CollectionFeature,
    Domain,
    Gene,
    Location,
    LocationPart,
    Module,
    Record,
)
from antismash_review.review import review_record

ROOT = Path(__file__).resolve().parent / "fixtures"


def _minimal_record(**kwargs: object) -> Record:
    defaults: dict[str, object] = {
        "record_id": "TEST.1",
        "name": "TEST",
        "description": "test record",
        "length": 100,
        "molecule_type": "DNA",
        "topology": "linear",
        "source_path": Path("/synthetic/test.gb"),
        "source_sha256": "0" * 64,
        "antismash_version": "8.0.4",
        "organism": "synthetic",
        "taxonomy": [],
    }
    defaults.update(kwargs)
    return Record(**defaults)  # type: ignore[arg-type]


def _simple_location(start: int, end: int, strand: int = 1) -> Location:
    part = LocationPart(start=start, end=end, strand=strand)
    return Location(
        start=start,
        end=end,
        strand=strand,
        parts=(part,),
        cross_origin=False,
        original=f"{start + 1}..{end}",
    )


def test_orphan_module_locus_diagnostic() -> None:
    """Module references a locus tag (G4) absent from the CDS set."""
    record = _minimal_record()
    record.genes.append(
        Gene(
            location=_simple_location(0, 30),
            locus_tag="G3",
            gene=None,
            product="synthetase",
            protein_id=None,
            translation=None,
            gene_kind="biosynthetic",
            gene_functions=[],
            ec_numbers=[],
            db_xrefs=[],
            notes=[],
            inference=[],
            region_numbers=[],
            candidate_cluster_numbers=[],
            protocluster_numbers=[],
            proto_core_numbers=[],
            qualifiers={},
        )
    )
    record.modules.append(
        Module(
            location=_simple_location(0, 60),
            domain_ids=[],
            locus_tags=["G3", "G4"],
            module_type="nrps",
            complete=True,
            starter=False,
            final=False,
            iterative=False,
            monomer_pairings=[],
            multi_cds=True,
        )
    )
    diagnostics = review_record(record)
    orphan = [d for d in diagnostics if d.code == "orphan_module_locus"]
    assert len(orphan) == 1
    assert "G4" in orphan[0].message


def test_no_orphan_when_all_locus_tags_present() -> None:
    """No orphan diagnostic when all module locus tags exist in the CDS set."""
    record = _minimal_record()
    record.genes.append(
        Gene(
            location=_simple_location(0, 30),
            locus_tag="G3",
            gene=None,
            product=None,
            protein_id=None,
            translation=None,
            gene_kind="biosynthetic",
            gene_functions=[],
            ec_numbers=[],
            db_xrefs=[],
            notes=[],
            inference=[],
            region_numbers=[],
            candidate_cluster_numbers=[],
            protocluster_numbers=[],
            proto_core_numbers=[],
            qualifiers={},
        )
    )
    record.modules.append(
        Module(
            location=_simple_location(0, 60),
            domain_ids=[],
            locus_tags=["G3"],
            module_type="nrps",
            complete=True,
            starter=False,
            final=False,
            iterative=False,
            monomer_pairings=[],
            multi_cds=False,
        )
    )
    diagnostics = review_record(record)
    assert not any(d.code == "orphan_module_locus" for d in diagnostics)


def test_missing_nrps_pks_architecture_diagnostic() -> None:
    """NRPS/PKS product advertised but no nrps_pks_domains domains present."""
    record = _minimal_record()
    record.regions.append(
        CollectionFeature(
            feature_type="region",
            number=1,
            location=_simple_location(0, 100),
            products=["NRPS"],
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
    diagnostics = review_record(record)
    missing = [d for d in diagnostics if d.code == "missing_nrps_pks_architecture"]
    assert len(missing) == 1
    assert "NRPS" in missing[0].message


def test_no_missing_architecture_when_domains_present() -> None:
    """No architecture diagnostic when nrps_pks_domains domains exist."""
    record = _minimal_record()
    record.regions.append(
        CollectionFeature(
            feature_type="region",
            number=1,
            location=_simple_location(0, 100),
            products=["NRPS"],
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
    record.domains.append(
        Domain(
            domain_id="d1",
            name="AMP-binding",
            subtypes=[],
            tool="nrps_pks_domains",
            locus_tag="G1",
            score=100.0,
            evalue=1e-50,
            protein_start=1,
            protein_end=100,
            specificity=[],
            location=_simple_location(10, 50),
            qualifiers={},
        )
    )
    diagnostics = review_record(record)
    assert not any(d.code == "missing_nrps_pks_architecture" for d in diagnostics)


def test_semantics_fixture_fires_orphan_module_diagnostic() -> None:
    """The semantics fixture's G4 locus tag should trigger orphan_module_locus."""
    record = parse_genbank(ROOT / "semantics.gb")[0]
    diagnostics = review_record(record)
    orphan = [d for d in diagnostics if d.code == "orphan_module_locus"]
    assert len(orphan) == 1
    assert "G4" in orphan[0].message


def test_partial_edge_fixture_fires_partial_cds() -> None:
    """The partial-edge fixture's <1..30 CDS should trigger partial_cds_at_edge."""
    record = parse_genbank(ROOT / "partial-edge.gb")[0]
    diagnostics = review_record(record)
    assert any(d.code == "partial_cds_at_edge" for d in diagnostics)
