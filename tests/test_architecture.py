from __future__ import annotations

from antismash_review.architecture import assess_architecture
from antismash_review.review import review_record
from tests.test_assemblyline import _domain, _module, _record


def test_canonical_t1pks_requires_ks_at_and_carrier_evidence() -> None:
    record = _record(
        [_module(100, 300, ["PKS1"], ["D1", "D2", "D3"], module_type="pks")],
        [
            _domain("D1", "PKS_KS", "PKS1"),
            _domain("D2", "PKS_AT", "PKS1"),
            _domain("D3", "PCP", "PKS1"),
        ],
        products=["T1PKS"],
    )

    assessment = assess_architecture(record)[0]

    assert assessment.status == "complete"
    assert assessment.score == 1.0
    assert assessment.missing_slots == ()


def test_t1pks_missing_at_is_partial_and_scoped_to_region() -> None:
    record = _record(
        [_module(100, 300, ["PKS1"], ["D1", "D2"], module_type="pks")],
        [_domain("D1", "PKS_KS", "PKS1"), _domain("D2", "PCP", "PKS1")],
        products=["T1PKS"],
    )

    assessment = assess_architecture(record)[0]

    assert assessment.scope == "region:1"
    assert assessment.status == "partial"
    assert assessment.missing_slots == ("AT",)
    assert assessment.score == 2 / 3


def test_trans_at_pks_exempts_missing_cis_at() -> None:
    record = _record(
        [_module(100, 300, ["PKS1"], ["D1", "D2"], module_type="pks")],
        [_domain("D1", "PKS_KS", "PKS1"), _domain("D2", "PCP", "PKS1")],
        products=["transAT-PKS"],
    )

    assessment = assess_architecture(record)[0]

    assert assessment.status == "complete"
    assert assessment.missing_slots == ()
    assert "cis-AT" in " ".join(assessment.exemptions)


def test_nrps_starter_and_internal_module_expectations_differ() -> None:
    starter = _record(
        [_module(100, 300, ["NRPS1"], ["D1", "D2"], starter=True)],
        [_domain("D1", "AMP-binding", "NRPS1"), _domain("D2", "PCP", "NRPS1")],
        products=["NRPS"],
    )
    internal = _record(
        [_module(100, 300, ["NRPS1"], ["D1", "D2", "D3"])],
        [
            _domain("D1", "Condensation", "NRPS1"),
            _domain("D2", "AMP-binding", "NRPS1"),
            _domain("D3", "PCP", "NRPS1"),
        ],
        products=["NRPS"],
    )

    starter_assessment = assess_architecture(starter)[0]
    internal_assessment = assess_architecture(internal)[0]

    assert starter_assessment.status == "complete"
    assert starter_assessment.expected_slots == ("A", "ACP/PCP")
    assert internal_assessment.status == "complete"
    assert internal_assessment.expected_slots == ("C", "A", "ACP/PCP")


def test_internal_nrps_without_condensation_is_partial() -> None:
    record = _record(
        [_module(100, 300, ["NRPS1"], ["D1", "D2"])],
        [_domain("D1", "AMP-binding", "NRPS1"), _domain("D2", "PCP", "NRPS1")],
        products=["NRPS"],
    )

    assessment = assess_architecture(record)[0]

    assert assessment.status == "partial"
    assert assessment.missing_slots == ("C",)


def test_incomplete_module_is_ambiguous_and_unsupported_label_is_not_applicable() -> None:
    incomplete = _record(
        [_module(100, 300, ["NRPS1"], ["D1", "D2"], complete=False)],
        [_domain("D1", "Condensation", "NRPS1"), _domain("D2", "AMP-binding", "NRPS1")],
        products=["NRPS"],
    )
    unsupported = _record(
        [_module(100, 300, ["PKS1"], ["D1"], module_type="pks")],
        [_domain("D1", "PKS_KS", "PKS1")],
        products=["NRPS-like"],
    )

    incomplete_assessment = assess_architecture(incomplete)[0]
    unsupported_assessment = assess_architecture(unsupported)[0]

    assert incomplete_assessment.status == "ambiguous"
    assert "antiSMASH marks" in " ".join(incomplete_assessment.caveats)
    assert unsupported_assessment.status == "not_applicable"


def test_edge_context_is_a_caveat_and_specific_warning_keeps_legacy_warning() -> None:
    record = _record(
        [_module(100, 300, ["PKS1"], ["D1"], module_type="pks")],
        [_domain("D1", "PKS_KS", "PKS1")],
        products=["T1PKS"],
    )
    record.regions[0].contig_edge = True

    assessment = assess_architecture(record)[0]
    codes = {diagnostic.code for diagnostic in review_record(record)}

    assert "assembly truncation" in " ".join(assessment.caveats)
    assert "architecture_core_domain_missing" in codes
    assert "missing_nrps_pks_architecture" not in codes
