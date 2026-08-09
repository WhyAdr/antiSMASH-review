"""Conservative domain-composition assessments scoped to regions and modules."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal

from .assemblyline import domains_for_module
from .locations import domains_in_collection, modules_in_collection
from .models import Domain, Record

ArchitectureStatus = Literal["complete", "partial", "ambiguous", "not_applicable"]


@dataclass(slots=True, frozen=True)
class DomainSlot:
    name: str
    aliases: tuple[str, ...]
    required: bool = True
    evidence_scope: str = "module"


@dataclass(slots=True, frozen=True)
class ArchitectureExpectation:
    name: str
    product_keys: tuple[str, ...]
    required_slots: tuple[DomainSlot, ...]
    optional_slots: tuple[DomainSlot, ...] = ()
    notes: tuple[str, ...] = ()


@dataclass(slots=True, frozen=True)
class ArchitectureAssessment:
    product: str
    scope: str
    region_number: int | None
    module_index: int | None
    status: ArchitectureStatus
    score: float | None
    expected_slots: tuple[str, ...]
    observed_slots: tuple[str, ...]
    missing_slots: tuple[str, ...]
    exemptions: tuple[str, ...]
    evidence_domains: tuple[str, ...]
    caveats: tuple[str, ...]


def _domain_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().casefold()).strip("_")


DOMAIN_FAMILIES: dict[str, tuple[str, ...]] = {
    "KS": ("ks", "pks_ks", "ketosynthase"),
    "AT": ("at", "pks_at", "acyltransferase"),
    "ACP/PCP": (
        "acp",
        "pcp",
        "pp_binding",
        "phosphopantetheine_binding",
        "carrier_protein",
    ),
    "A": ("a", "amp_binding", "adenylation", "adenylation_domain"),
    "C": ("c", "condensation", "condensation_domain"),
    "TE": ("te", "thioesterase", "thioesterase_like", "thioesterase_domain"),
}


_FAMILY_KEYS = {
    family: frozenset(_domain_key(alias) for alias in aliases)
    for family, aliases in DOMAIN_FAMILIES.items()
}


def domain_families(domains: list[Domain]) -> tuple[str, ...]:
    """Return recognized aSDomain families in stable registry order."""

    observed = {
        _domain_key(value)
        for domain in domains
        for value in ((domain.name,) if domain.name else ()) + tuple(domain.subtypes)
    }
    return tuple(
        family for family in DOMAIN_FAMILIES if observed.intersection(_FAMILY_KEYS[family])
    )


def _has_family(domains: list[Domain], family: str) -> bool:
    return family in domain_families(domains)


def _score_and_status(
    expected: tuple[DomainSlot, ...],
    domains: list[Domain],
    *,
    ambiguous: bool = False,
) -> tuple[ArchitectureStatus, float | None, tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    expected_names = tuple(slot.name for slot in expected if slot.required)
    missing = tuple(
        slot.name for slot in expected if slot.required and not _has_family(domains, slot.name)
    )
    satisfied = len(expected_names) - len(missing)
    score = satisfied / len(expected_names) if expected_names else None
    if ambiguous:
        status: ArchitectureStatus = "ambiguous"
    elif missing:
        status = "partial"
    else:
        status = "complete"
    return status, score, expected_names, tuple(missing), domain_families(domains)


def _assessment(
    *,
    product: str,
    scope: str,
    region_number: int | None,
    module_index: int | None,
    expected: tuple[DomainSlot, ...],
    domains: list[Domain],
    exemptions: tuple[str, ...] = (),
    caveats: tuple[str, ...] = (),
    ambiguous: bool = False,
) -> ArchitectureAssessment:
    status, score, expected_names, missing, observed = _score_and_status(
        expected,
        domains,
        ambiguous=ambiguous,
    )
    return ArchitectureAssessment(
        product=product,
        scope=scope,
        region_number=region_number,
        module_index=module_index,
        status=status,
        score=score,
        expected_slots=expected_names,
        observed_slots=observed,
        missing_slots=missing,
        exemptions=exemptions,
        evidence_domains=tuple(domain.name for domain in domains if domain.name is not None),
        caveats=(
            "score measures expected parsed-domain coverage, not biosynthetic activity",
            *caveats,
        ),
    )


_T1PKS = ArchitectureExpectation(
    name="canonical T1PKS",
    product_keys=("t1pks",),
    required_slots=(
        DomainSlot("KS", DOMAIN_FAMILIES["KS"], evidence_scope="region"),
        DomainSlot("AT", DOMAIN_FAMILIES["AT"], evidence_scope="region"),
        DomainSlot("ACP/PCP", DOMAIN_FAMILIES["ACP/PCP"], evidence_scope="region"),
    ),
)
_TRANS_AT_PKS = ArchitectureExpectation(
    name="trans-AT PKS",
    product_keys=("transat_pks", "trans_at_pks"),
    required_slots=(
        DomainSlot("KS", DOMAIN_FAMILIES["KS"], evidence_scope="region"),
        DomainSlot("ACP/PCP", DOMAIN_FAMILIES["ACP/PCP"], evidence_scope="region"),
    ),
    notes=("cis-AT absence is expected for trans-AT product labels",),
)


def _product_key(product: str) -> str:
    return _domain_key(product)


def _region_edge_caveat(record: Record, region_index: int) -> str | None:
    region = record.regions[region_index]
    if (
        region.contig_edge is True
        or region.location.start == 0
        or region.location.end == record.length
    ):
        return "missing core-domain evidence may reflect record or assembly truncation"
    return None


def _region_assessment(
    record: Record,
    region_index: int,
    product: str,
    expectation: ArchitectureExpectation,
) -> ArchitectureAssessment:
    region = record.regions[region_index]
    domains = domains_in_collection(record, region)
    caveats = tuple(
        item for item in (_region_edge_caveat(record, region_index),) if item is not None
    )
    exemptions = expectation.notes if expectation.name == "trans-AT PKS" else ()
    return _assessment(
        product=product,
        scope=f"region:{region.number if region.number is not None else region_index + 1}",
        region_number=region.number,
        module_index=None,
        expected=expectation.required_slots,
        domains=domains,
        exemptions=exemptions,
        caveats=caveats,
    )


def _nrps_assessments(
    record: Record, region_index: int, product: str
) -> list[ArchitectureAssessment]:
    region = record.regions[region_index]
    modules = modules_in_collection(record, region)
    if not modules:
        return [
            ArchitectureAssessment(
                product=product,
                scope=f"region:{region.number if region.number is not None else region_index + 1}",
                region_number=region.number,
                module_index=None,
                status="not_applicable",
                score=None,
                expected_slots=(),
                observed_slots=(),
                missing_slots=(),
                exemptions=(),
                evidence_domains=(),
                caveats=("no explicit antiSMASH modules were available for an NRPS assessment",),
            )
        ]

    assessments: list[ArchitectureAssessment] = []
    for module_index, module in enumerate(modules, start=1):
        if module.module_type is not None and module.module_type.casefold() != "nrps":
            continue
        domains = domains_for_module(record, module)
        expected = [
            DomainSlot("A", DOMAIN_FAMILIES["A"]),
            DomainSlot("ACP/PCP", DOMAIN_FAMILIES["ACP/PCP"]),
        ]
        if not module.starter:
            expected.insert(0, DomainSlot("C", DOMAIN_FAMILIES["C"]))
        caveats: list[str] = []
        edge = _region_edge_caveat(record, region_index)
        if edge is not None:
            caveats.append(edge)
        if not module.complete:
            caveats.append("antiSMASH marks this module incomplete")
        if module.final:
            caveats.append(
                "antiSMASH marks this module final; termination chemistry is not inferred"
            )
        if module.iterative:
            caveats.append("antiSMASH marks this module iterative; residue count is not inferred")
        scope = (
            f"region:{region.number if region.number is not None else region_index + 1}"
            f"/module:{module_index}"
        )
        assessments.append(
            _assessment(
                product=product,
                scope=scope,
                region_number=region.number,
                module_index=module_index,
                expected=tuple(expected),
                domains=domains,
                caveats=tuple(caveats),
                ambiguous=not module.complete,
            )
        )
    return assessments


def assess_architecture(record: Record) -> list[ArchitectureAssessment]:
    """Assess only conservative, source-auditable product/domain expectations."""

    assessments: list[ArchitectureAssessment] = []
    for region_index, region in enumerate(record.regions):
        for product in region.products:
            key = _product_key(product)
            if key in _T1PKS.product_keys:
                assessments.append(_region_assessment(record, region_index, product, _T1PKS))
            elif key in _TRANS_AT_PKS.product_keys:
                assessments.append(_region_assessment(record, region_index, product, _TRANS_AT_PKS))
            elif key == "nrps":
                assessments.extend(_nrps_assessments(record, region_index, product))
            else:
                scope = f"region:{region.number if region.number is not None else region_index + 1}"
                assessments.append(
                    ArchitectureAssessment(
                        product=product,
                        scope=scope,
                        region_number=region.number,
                        module_index=None,
                        status="not_applicable",
                        score=None,
                        expected_slots=(),
                        observed_slots=(),
                        missing_slots=(),
                        exemptions=(),
                        evidence_domains=(),
                        caveats=(
                            "no conservative Phase 3 domain expectation is registered "
                            "for this product class",
                        ),
                    )
                )
    return assessments
