from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import pytest

from antismash_review.assemblyline import (
    domains_for_module,
    predict_assembly_lines,
)
from antismash_review.exporters.assemblyline_json import dumps_assembly_lines
from antismash_review.exporters.assemblyline_markdown import render_assemblyline_markdown
from antismash_review.exporters.assemblyline_table import render_assemblyline_tsv
from antismash_review.models import (
    CollectionFeature,
    Domain,
    Location,
    LocationPart,
    Module,
    Record,
)


def _location(start: int, end: int, strand: int | None = 1) -> Location:
    return Location(
        start=start,
        end=end,
        strand=strand,
        parts=(LocationPart(start, end, strand),),
        cross_origin=False,
        original=f"{start}..{end}",
    )


def _domain(
    domain_id: str,
    name: str,
    locus_tag: str,
    *,
    protein_start: int | None = None,
    protein_end: int | None = None,
    specificity: list[str] | None = None,
) -> Domain:
    return Domain(
        domain_id=domain_id,
        name=name,
        subtypes=[],
        tool="nrps_pks_domains",
        locus_tag=locus_tag,
        score=None,
        evalue=None,
        protein_start=protein_start,
        protein_end=protein_end,
        specificity=specificity or [],
        location=_location(0, 100),
        qualifiers={},
    )


def _module(
    start: int,
    end: int,
    locus_tags: list[str],
    domain_ids: list[str],
    *,
    strand: int | None = 1,
    pairings: list[str] | None = None,
    module_type: str | None = "nrps",
    complete: bool = True,
    starter: bool = False,
    final: bool = False,
    iterative: bool = False,
) -> Module:
    return Module(
        location=_location(start, end, strand),
        domain_ids=domain_ids,
        locus_tags=locus_tags,
        module_type=module_type,
        complete=complete,
        starter=starter,
        final=final,
        iterative=iterative,
        monomer_pairings=pairings or [],
        multi_cds=len(locus_tags) > 1,
    )


def _record(
    modules: list[Module],
    domains: list[Domain],
    *,
    products: list[str] | None = None,
) -> Record:
    return Record(
        record_id="ASSEMBLY.1",
        name="ASSEMBLY.1",
        description="synthetic assembly-line record",
        length=1000,
        molecule_type="DNA",
        topology="linear",
        source_path=Path("assembly.gbk"),
        source_sha256="",
        antismash_version="8.0.4",
        organism=None,
        taxonomy=[],
        regions=[
            CollectionFeature(
                feature_type="region",
                number=1,
                location=_location(0, 1000),
                products=products or ["NRPS"],
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
        modules=modules,
        domains=domains,
    )


def test_forward_modules_use_protein_order_and_retain_all_pairings() -> None:
    domains = [
        _domain("D1", "AMP-binding", "CDS1", protein_start=10, protein_end=30),
        _domain("D2", "AMP-binding", "CDS1", protein_start=110, protein_end=130),
        _domain("D3", "AMP-binding", "CDS1", protein_start=210, protein_end=230),
    ]
    record = _record(
        [
            _module(300, 400, ["CDS1"], ["D3"], pairings=["Thr -> Thr"], final=True),
            _module(100, 200, ["CDS1"], ["D1"], pairings=["Ser -> Ser"], starter=True),
            _module(200, 300, ["CDS1"], ["D2"], pairings=["Leu -> Leu"]),
        ],
        domains,
    )

    prediction = predict_assembly_lines(record)[0]

    assert prediction.ordering_basis == "protein-domain-order"
    assert [call.monomer for call in prediction.chain] == ["Ser", "Leu", "Thr"]
    assert prediction.modules[0].starter is True
    assert prediction.modules[-1].final is True
    assert "mass" in asdict(prediction)
    assert prediction.mass is not None
    assert prediction.mass.linear_core_mass_da is None


def test_reverse_strand_modules_use_reverse_genomic_order_without_protein_coordinates() -> None:
    record = _record(
        [
            _module(100, 200, ["CDS-R"], ["D1"], strand=-1, pairings=["Ser -> Ser"]),
            _module(300, 400, ["CDS-R"], ["D2"], strand=-1, pairings=["Leu -> Leu"]),
        ],
        [_domain("D1", "AMP-binding", "CDS-R"), _domain("D2", "AMP-binding", "CDS-R")],
    )

    prediction = predict_assembly_lines(record)[0]

    assert prediction.ordering_basis == "strand-aware-nucleotide-order"
    assert [call.monomer for call in prediction.chain] == ["Leu", "Ser"]
    assert prediction.strand == -1


def test_multi_cds_module_is_not_split_or_merged_by_proximity() -> None:
    record = _record(
        [
            _module(
                100,
                200,
                ["CDS-A", "CDS-B"],
                ["D1", "D2"],
                pairings=["X -> X", "Ser -> Ser"],
            ),
            _module(220, 320, ["CDS-A"], ["D3"], pairings=["Leu -> Leu"]),
        ],
        [
            _domain("D1", "Condensation", "CDS-A"),
            _domain("D2", "AMP-binding", "CDS-B"),
            _domain("D3", "AMP-binding", "CDS-A"),
        ],
    )

    predictions = predict_assembly_lines(record)
    multi = next(prediction for prediction in predictions if prediction.modules[0].multi_cds)

    assert len(predictions) == 2
    assert len(multi.modules) == 1
    assert [call.monomer for call in multi.chain] == ["X", "Ser"]
    assert "multi-CDS membership" in " ".join(multi.caveats)
    assert multi.mass is not None
    assert multi.mass.linear_core_mass_da is None
    assert "X" in multi.mass.unresolved_monomers


def test_distinct_cds_chains_in_one_region_are_not_collapsed() -> None:
    record = _record(
        [
            _module(100, 200, ["CDS-A"], ["D1"], pairings=["Ser -> Ser"]),
            _module(210, 310, ["CDS-B"], ["D2"], pairings=["Leu -> Leu"]),
        ],
        [_domain("D1", "AMP-binding", "CDS-A"), _domain("D2", "AMP-binding", "CDS-B")],
    )

    predictions = predict_assembly_lines(record)

    assert len(predictions) == 2
    assert {prediction.chain[0].monomer for prediction in predictions} == {"Ser", "Leu"}


def test_specificity_consensus_is_a_low_confidence_fallback() -> None:
    record = _record(
        [_module(100, 200, ["CDS1"], ["D1"])],
        [
            _domain(
                "D1",
                "AMP-binding",
                "CDS1",
                specificity=["substrate consensus: Ser", "Minowa: Ser"],
            )
        ],
    )

    call = predict_assembly_lines(record)[0].chain[0]

    assert call.monomer == "Ser"
    assert call.source == "domain_specificity"
    assert call.confidence == "low"
    assert any("Minowa: Ser" in note for note in call.notes)


def test_conflicting_specificity_outputs_are_retained_as_ambiguous() -> None:
    record = _record(
        [_module(100, 200, ["CDS1"], ["D1"])],
        [
            _domain(
                "D1",
                "AMP-binding",
                "CDS1",
                specificity=["substrate consensus: Ser", "substrate consensus: Leu"],
            )
        ],
    )

    prediction = predict_assembly_lines(record)[0]

    assert [call.monomer for call in prediction.chain] == ["Ser", "Leu"]
    assert all(call.confidence == "unresolved" for call in prediction.chain)
    assert "conflicting substrate specificity" in " ".join(prediction.modules[0].warnings)


def test_nonproteinogenic_pairing_is_preserved_without_amino_acid_normalization() -> None:
    record = _record(
        [_module(100, 200, ["CDS1"], ["D1"], pairings=["mal -> ccmal"])],
        [_domain("D1", "PKS_KS", "CDS1")],
    )

    call = predict_assembly_lines(record)[0].chain[0]

    assert call.substrate == "mal"
    assert call.monomer == "ccmal"
    assert call.confidence == "high"


def test_fully_resolved_peptide_has_separate_linear_and_cyclic_candidates() -> None:
    record = _record(
        [
            _module(100, 200, ["CDS1"], ["D1"], pairings=["Ser -> Ser"]),
            _module(210, 310, ["CDS1"], ["D2"], pairings=["Leu -> Leu"]),
        ],
        [
            _domain("D1", "AMP-binding", "CDS1", protein_start=10, protein_end=20),
            _domain("D2", "AMP-binding", "CDS1", protein_start=30, protein_end=40),
        ],
    )

    mass = predict_assembly_lines(record)[0].mass
    assert mass is not None

    water = 2 * 1.00782503223 + 15.99491461957
    # Free Ser: C3H7NO3, Free Leu: C6H13NO2
    ser = 3 * 12.0 + 7 * 1.00782503223 + 14.00307400443 + 3 * 15.99491461957
    leu = 6 * 12.0 + 13 * 1.00782503223 + 14.00307400443 + 2 * 15.99491461957
    expected_linear = ser + leu - water
    assert mass.linear_core_mass_da == pytest.approx(expected_linear, abs=1e-9)
    assert mass.head_to_tail_cyclic_candidate_mass_da == pytest.approx(
        expected_linear - water,
        abs=1e-9,
    )
    assert mass.resolved_monomers == 2
    assert mass.total_monomers == 2
    assert mass.coverage_fraction == 1.0
    assert mass.topology_assumption == "unknown"


def test_starter_module_blocks_mass_when_tail_chemistry_is_unresolved() -> None:
    record = _record(
        [_module(100, 200, ["CDS1"], ["D1"], pairings=["Ser -> Ser"], starter=True)],
        [_domain("D1", "AMP-binding", "CDS1")],
    )

    mass = predict_assembly_lines(record)[0].mass

    assert mass is not None
    assert mass.linear_core_mass_da is None
    assert "starter chemistry" in mass.chemistry_scope


def test_hybrid_pks_segment_blocks_peptide_only_mass() -> None:
    record = _record(
        [
            _module(100, 200, ["CDS1"], ["D1"], pairings=["Ser -> Ser"]),
            _module(
                210,
                310,
                ["CDS1"],
                ["D2"],
                pairings=["mal -> ccmal"],
                module_type="pks",
            ),
        ],
        [
            _domain("D1", "AMP-binding", "CDS1", protein_start=10, protein_end=20),
            _domain("D2", "PKS_KS", "CDS1", protein_start=30, protein_end=40),
        ],
    )

    prediction = predict_assembly_lines(record)[0]

    assert prediction.mass is not None
    assert prediction.mass.linear_core_mass_da is None
    assert "non-NRPS" in prediction.mass.chemistry_scope
    assert "PKS evidence" in " ".join(prediction.caveats)


def test_release_domain_and_flags_are_reported_without_release_mode_claim() -> None:
    record = _record(
        [
            _module(
                100,
                200,
                ["CDS1"],
                ["D1", "D2"],
                pairings=["Ser -> Ser"],
                final=True,
                iterative=True,
            )
        ],
        [
            _domain("D1", "AMP-binding", "CDS1"),
            _domain("D2", "Thioesterase", "CDS1"),
        ],
    )

    prediction = predict_assembly_lines(record)[0]

    assert prediction.modules[0].release_domains == ("Thioesterase",)
    assert prediction.modules[0].final is True
    assert prediction.modules[0].iterative is True
    assert "release mode" in " ".join(prediction.caveats)
    assert "cyclization" in " ".join(prediction.caveats)


def test_module_domain_resolution_does_not_guess_duplicate_ids() -> None:
    module = _module(100, 200, ["CDS1"], ["D1", "MISSING"])
    record = _record(
        [module],
        [_domain("D1", "AMP-binding", "CDS1"), _domain("D1", "PCP", "CDS1")],
    )

    assert domains_for_module(record, module) == []


def test_assemblyline_exports_are_deterministic_and_versioned() -> None:
    record = _record(
        [_module(100, 200, ["CDS1"], ["D1"], pairings=["Ser -> Ser"])],
        [_domain("D1", "AMP-binding", "CDS1")],
    )

    json_output = dumps_assembly_lines([record])
    table_output = render_assemblyline_tsv([record])
    markdown_output = render_assemblyline_markdown([record])

    assert '"schema_name": "antismash-review-assemblyline"' in json_output
    assert "record_id\tregion_number\tassembly_line" in table_output
    assert "Ser" in markdown_output
    assert json_output == dumps_assembly_lines([record])
    assert table_output == render_assemblyline_tsv([record])
