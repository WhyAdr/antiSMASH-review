from __future__ import annotations

from pathlib import Path

from antismash_review.fingerprints import domain_counter, product_counter
from antismash_review.models import CollectionFeature, Location, LocationPart, Record


def _record() -> Record:
    location = Location(0, 100, 1, (LocationPart(0, 100, 1),), False, "0..100")
    return Record(
        record_id="FP.1",
        name="FP.1",
        description="fingerprint fixture",
        length=100,
        molecule_type="DNA",
        topology="linear",
        source_path=Path("fp.gbk"),
        source_sha256="",
        antismash_version=None,
        organism=None,
        taxonomy=[],
        regions=[
            CollectionFeature(
                feature_type="region",
                number=1,
                location=location,
                products=[" NRPS ", "NRPS"],
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
                contig_edge=False,
                qualifiers={},
            )
        ],
    )


def test_fingerprints_keep_raw_tokens_separate_from_normalized_keys() -> None:
    record = _record()

    assert product_counter([record]) == {" NRPS ": 1, "NRPS": 1}
    assert product_counter([record], normalized=True) == {"nrps": 2}
    assert domain_counter([record]) == {}
