"""Evidence-only interpretation of antiSMASH NRPS/PKS modules.

This module deliberately stops at antiSMASH-derived module and monomer evidence.  It
does not calculate chemistry, infer a final metabolite, or join separate CDS-local
chains using genomic proximity.
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from .chemistry import WATER_MASS, free_monomer_mass
from .models import Domain, Module, Record

CallSource = Literal["module_pairing", "domain_specificity", "unknown"]
CallConfidence = Literal["high", "medium", "low", "unresolved"]
OrderingConfidence = Literal["high", "medium", "low"]

_PAIRING_RE = re.compile(r"^\s*(?P<substrate>.*?)\s*->\s*(?P<monomer>.*?)\s*$")
_CONSENSUS_RE = re.compile(r"^\s*substrate\s+consensus\s*:\s*(?P<value>.+?)\s*$", re.IGNORECASE)
_UNRESOLVED_MONOMERS = frozenset({"", "?", "x", "unknown", "unresolved"})
_ADENYLATION_NAMES = frozenset(
    {
        "a",
        "amp-binding",
        "amp_binding",
        "adenylation",
        "adenylation_domain",
    }
)
_RELEASE_DOMAIN_ALIASES = frozenset(
    {
        "te",
        "thioesterase",
        "thioesterase-like",
        "thioesterase_like",
        "thioesterase domain",
        "thioesterase_domain",
    }
)


@dataclass(slots=True, frozen=True)
class MonomerCall:
    """One antiSMASH-derived substrate/monomer call, including unresolved calls."""

    substrate: str | None
    monomer: str | None
    display: str
    source: CallSource
    confidence: CallConfidence
    notes: tuple[str, ...] = ()


@dataclass(slots=True, frozen=True)
class ModulePrediction:
    """Resolved evidence and calls for one antiSMASH module."""

    index: int
    module_type: str | None
    locus_tags: tuple[str, ...]
    domain_ids: tuple[str, ...]
    domain_names: tuple[str, ...]
    complete: bool
    starter: bool
    final: bool
    iterative: bool
    multi_cds: bool
    monomer_calls: tuple[MonomerCall, ...]
    release_domains: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(slots=True, frozen=True)
class MassEstimate:
    """Conservative candidate masses for a fully resolved amino-acid-like core."""

    linear_core_mass_da: float | None
    head_to_tail_cyclic_candidate_mass_da: float | None
    resolved_monomers: int
    total_monomers: int
    coverage_fraction: float
    unresolved_monomers: tuple[str, ...]
    topology_assumption: Literal["linear", "cyclic", "unknown"]
    chemistry_scope: str


@dataclass(slots=True, frozen=True)
class AssemblyLinePrediction:
    """A deterministic local assembly-line hypothesis and optional core-mass candidates."""

    record_id: str
    region_number: int | None
    strand: int | None
    modules: tuple[ModulePrediction, ...]
    chain: tuple[MonomerCall, ...]
    ordering_basis: str
    ordering_confidence: OrderingConfidence
    caveats: tuple[str, ...]
    mass: MassEstimate | None = None


def _normalized_domain_name(name: str) -> str:
    return name.strip().casefold().replace(" ", "_")


def _is_adenylation_domain(domain: Domain) -> bool:
    if domain.name is None:
        return False
    return _normalized_domain_name(domain.name) in _ADENYLATION_NAMES


def _is_unresolved(value: str | None) -> bool:
    return value is None or value.strip().casefold() in _UNRESOLVED_MONOMERS


def _domain_is_release(domain: Domain) -> bool:
    if domain.name is None:
        return False
    return _normalized_domain_name(domain.name) in _RELEASE_DOMAIN_ALIASES


def domains_for_module(record: Record, module: Module) -> list[Domain]:
    """Resolve a module's domain IDs without guessing duplicated domain instances."""

    by_id: dict[str, list[Domain]] = defaultdict(list)
    for domain in record.domains:
        if domain.domain_id is not None:
            by_id[domain.domain_id].append(domain)

    resolved: list[Domain] = []
    for domain_id in module.domain_ids:
        matches = by_id.get(domain_id, [])
        if len(matches) == 1:
            resolved.append(matches[0])
    return resolved


def _pairing_call(raw: str) -> MonomerCall:
    match = _PAIRING_RE.match(raw)
    if match is None:
        return MonomerCall(
            substrate=None,
            monomer=None,
            display="?",
            source="module_pairing",
            confidence="unresolved",
            notes=(
                f"raw module pairing retained: {raw}",
                "malformed pairing; expected substrate -> monomer",
            ),
        )

    substrate = match.group("substrate").strip() or None
    monomer = match.group("monomer").strip() or None
    confidence: CallConfidence = "unresolved" if _is_unresolved(monomer) else "high"
    notes = [f"raw module pairing retained: {raw}"]
    if confidence == "unresolved":
        notes.append("monomer is unresolved in the antiSMASH call")
    return MonomerCall(
        substrate=substrate,
        monomer=monomer,
        display=monomer or "?",
        source="module_pairing",
        confidence=confidence,
        notes=tuple(notes),
    )


def _specificity_calls(
    domains: Iterable[Domain],
) -> tuple[tuple[MonomerCall, ...], tuple[str, ...]]:
    adenylation_domains = [domain for domain in domains if _is_adenylation_domain(domain)]
    raw_specificity = [value for domain in adenylation_domains for value in domain.specificity]
    candidates: list[tuple[str, str]] = []
    for raw in raw_specificity:
        match = _CONSENSUS_RE.match(raw)
        if match:
            value = match.group("value").strip()
            if value and all(value.casefold() != existing.casefold() for existing, _ in candidates):
                candidates.append((value, raw))

    if len(candidates) == 1:
        value, _raw = candidates[0]
        confidence: CallConfidence = "unresolved" if _is_unresolved(value) else "low"
        notes = tuple(f"raw specificity retained: {item}" for item in raw_specificity)
        if confidence == "unresolved":
            notes += ("specificity consensus is unresolved",)
        return (
            MonomerCall(
                substrate=None,
                monomer=value,
                display=value,
                source="domain_specificity",
                confidence=confidence,
                notes=notes,
            ),
        ), ()

    if len(candidates) > 1:
        calls = tuple(
            MonomerCall(
                substrate=None,
                monomer=value,
                display=value,
                source="domain_specificity",
                confidence="unresolved",
                notes=(
                    *(f"raw specificity retained: {item}" for item in raw_specificity),
                    "conflicting specificity outputs retained; no single call selected",
                ),
            )
            for value, _raw in candidates
        )
        return calls, ("conflicting substrate specificity outputs retained",)

    if raw_specificity:
        return (
            MonomerCall(
                substrate=None,
                monomer=None,
                display="?",
                source="domain_specificity",
                confidence="unresolved",
                notes=(
                    *(f"raw specificity retained: {item}" for item in raw_specificity),
                    "no supported substrate consensus form found",
                ),
            ),
        ), ("specificity evidence was retained but not interpreted",)

    return (
        MonomerCall(
            substrate=None,
            monomer=None,
            display="?",
            source="unknown",
            confidence="unresolved",
            notes=("no module pairing or supported specificity evidence",),
        ),
    ), ()


def _calls_for_module(
    module: Module,
    domains: list[Domain],
) -> tuple[tuple[MonomerCall, ...], list[str]]:
    warnings: list[str] = []
    if module.monomer_pairings:
        calls = tuple(_pairing_call(raw) for raw in module.monomer_pairings)
        if len(calls) > 1:
            warnings.append("multiple module pairing calls retained; no single call selected")
        if any(call.confidence == "unresolved" for call in calls):
            warnings.append("module pairing contains unresolved monomer evidence")
        return calls, warnings

    calls, specificity_warnings = _specificity_calls(domains)
    warnings.extend(specificity_warnings)
    return calls, warnings


def _module_prediction(module: Module, index: int, domains: list[Domain]) -> ModulePrediction:
    calls, warnings = _calls_for_module(module, domains)
    resolved_ids = {domain.domain_id for domain in domains if domain.domain_id is not None}
    missing_ids = [domain_id for domain_id in module.domain_ids if domain_id not in resolved_ids]
    if missing_ids:
        warnings.append(f"unresolved domain IDs retained: {', '.join(missing_ids)}")
    if not module.complete:
        warnings.append("antiSMASH module is marked incomplete")
    if module.multi_cds:
        warnings.append("multi-CDS module retained as one antiSMASH module")

    release_domains = tuple(
        domain.name for domain in domains if domain.name is not None and _domain_is_release(domain)
    )
    return ModulePrediction(
        index=index,
        module_type=module.module_type,
        locus_tags=tuple(module.locus_tags),
        domain_ids=tuple(module.domain_ids),
        domain_names=tuple(domain.name for domain in domains if domain.name is not None),
        complete=module.complete,
        starter=module.starter,
        final=module.final,
        iterative=module.iterative,
        multi_cds=module.multi_cds,
        monomer_calls=calls,
        release_domains=release_domains,
        warnings=tuple(warnings),
    )


def _region_indices(record: Record, module: Module) -> tuple[int, ...]:
    return tuple(
        index
        for index, region in enumerate(record.regions)
        if any(
            module_part.start < region_part.end and region_part.start < module_part.end
            for module_part in module.location.parts
            for region_part in region.location.parts
        )
    )


def _protein_coordinates(record: Record, module: Module) -> list[tuple[int, int]]:
    return [
        (domain.protein_start, domain.protein_end)
        for domain in domains_for_module(record, module)
        if domain.protein_start is not None and domain.protein_end is not None
    ]


def _ordered_modules(
    record: Record,
    indexed_modules: list[tuple[int, Module]],
    *,
    multi_cds_group: bool,
) -> tuple[list[tuple[int, Module]], str, OrderingConfidence]:
    if multi_cds_group:
        return indexed_modules, "antiSMASH-module-membership", "high"

    if all(_protein_coordinates(record, module) for _, module in indexed_modules):
        ordered = sorted(
            indexed_modules,
            key=lambda item: (
                min(start for start, _end in _protein_coordinates(record, item[1])),
                max(end for _start, end in _protein_coordinates(record, item[1])),
                item[0],
            ),
        )
        return ordered, "protein-domain-order", "high"

    strands = {module.location.strand for _, module in indexed_modules}
    if len(strands) == 1 and next(iter(strands)) in {-1, 1}:
        strand = next(iter(strands))
        if strand == -1:

            def key(item: tuple[int, Module]) -> tuple[int, int, int]:
                return -item[1].location.end, -item[1].location.start, item[0]
        else:

            def key(item: tuple[int, Module]) -> tuple[int, int, int]:
                return item[1].location.start, item[1].location.end, item[0]

        return sorted(indexed_modules, key=key), "strand-aware-nucleotide-order", "high"

    ordered = sorted(
        indexed_modules,
        key=lambda item: (item[1].location.start, item[1].location.end, item[0]),
    )
    return ordered, "genomic-coordinate-order", "low"


def _group_sort_key(
    group: tuple[tuple[str, tuple[int, ...], str], list[tuple[int, Module]]],
) -> tuple[int, int, str]:
    key, modules = group
    first = min(module.location.start for _index, module in modules)
    return first, min(key[1], default=-1), key[2]


def _region_number(record: Record, modules: Iterable[Module]) -> int | None:
    numbers = {
        record.regions[index].number
        for module in modules
        for index in _region_indices(record, module)
        if record.regions[index].number is not None
    }
    return next(iter(numbers)) if len(numbers) == 1 else None


def _strand(modules: Iterable[Module]) -> int | None:
    strands = {module.location.strand for module in modules}
    return next(iter(strands)) if len(strands) == 1 else None


def _prediction_caveats(
    modules: list[Module],
    module_predictions: list[ModulePrediction],
    *,
    ordering_basis: str,
    ordering_confidence: OrderingConfidence,
) -> tuple[str, ...]:
    caveats: list[str] = []
    if ordering_confidence != "high":
        caveats.append("module ordering is lower-confidence genomic evidence")
    if ordering_basis == "antiSMASH-module-membership":
        caveats.append("separate CDS-local order was not inferred for this multi-CDS module")
    if any(module.multi_cds for module in modules):
        caveats.append(
            "multi-CDS membership is retained from antiSMASH without a cross-CDS heuristic"
        )
    if any(module.iterative for module in modules):
        caveats.append("iterative module flags do not imply a one-module-one-incorporation count")
    if any(pred.release_domains for pred in module_predictions):
        caveats.append(
            "release-domain evidence is present; release mode is unknown; "
            "hydrolysis versus cyclization is unresolved"
        )
    if any(
        prediction.module_type is not None and prediction.module_type.casefold() == "pks"
        for prediction in module_predictions
    ):
        caveats.append(
            "PKS evidence is reported as a module call, not converted into peptide chemistry"
        )
    if any(
        call.confidence == "unresolved"
        for prediction in module_predictions
        for call in prediction.monomer_calls
    ):
        caveats.append("one or more monomer calls remain unresolved or ambiguous")
    return tuple(caveats)


def _estimate_core_mass(
    modules: tuple[ModulePrediction, ...],
    chain: tuple[MonomerCall, ...],
) -> MassEstimate:
    total = len(chain)
    resolved = 0
    masses: list[float] = []
    unresolved: list[str] = []
    reasons: list[str] = []

    for call in chain:
        mass = free_monomer_mass(call.monomer)
        if mass is None or call.confidence != "high":
            unresolved.append(call.monomer or call.display or "?")
            if call.confidence != "high":
                reasons.append("one or more calls are not high-confidence module pairings")
        else:
            resolved += 1
            masses.append(mass)

    for module in modules:
        module_type = module.module_type.casefold() if module.module_type else None
        if module_type != "nrps":
            reasons.append("non-NRPS or unknown module chemistry is present")
        if module.iterative:
            reasons.append("iterative module incorporation count is unresolved")
        if module.starter:
            reasons.append("starter chemistry may include an unresolved acyl or other tail")

    coverage = resolved / total if total else 0.0
    unique_reasons = tuple(dict.fromkeys(reasons))
    complete = total > 0 and resolved == total and not unique_reasons
    linear: float | None = None
    cyclic: float | None = None
    if complete:
        water_mass = WATER_MASS
        linear = sum(masses) - (total - 1) * water_mass
        cyclic = linear - water_mass

    if complete:
        scope = (
            "fully resolved proteinogenic amino-acid-like NRPS core; candidates exclude "
            "tailoring chemistry and observed metabolite identity"
        )
    else:
        scope = "full-core candidate unavailable: " + "; ".join(
            unique_reasons or ("unresolved chemistry",)
        )

    return MassEstimate(
        linear_core_mass_da=linear,
        head_to_tail_cyclic_candidate_mass_da=cyclic,
        resolved_monomers=resolved,
        total_monomers=total,
        coverage_fraction=coverage,
        unresolved_monomers=tuple(dict.fromkeys(unresolved)),
        topology_assumption="unknown",
        chemistry_scope=scope,
    )


def _mass_caveats(mass: MassEstimate) -> tuple[str, ...]:
    caveats: list[str] = []
    if mass.unresolved_monomers:
        caveats.append(
            "full-core mass is unavailable for unresolved or unmodeled monomers: "
            + ", ".join(mass.unresolved_monomers)
        )
    if mass.linear_core_mass_da is not None:
        caveats.append(
            "linear and cyclic values are modeled core candidates, not the measured "
            "final metabolite mass"
        )
    if mass.linear_core_mass_da is None and not mass.unresolved_monomers:
        caveats.append(mass.chemistry_scope)
    return tuple(caveats)


def predict_assembly_lines(record: Record) -> list[AssemblyLinePrediction]:
    """Build deterministic local assembly-line predictions from parsed antiSMASH modules."""

    groups: dict[tuple[str, tuple[int, ...], str], list[tuple[int, Module]]] = defaultdict(list)
    for index, module in enumerate(record.modules):
        scope = _region_indices(record, module) or (-1,)
        if len(module.locus_tags) == 1:
            group_name = f"cds:{module.locus_tags[0]}"
            group_type = "cds"
        elif module.multi_cds or len(module.locus_tags) > 1:
            group_name = f"module:{index}"
            group_type = "multi-cds"
        else:
            group_name = f"module:{index}"
            group_type = "unassigned"
        groups[(group_type, scope, group_name)].append((index, module))

    predictions: list[AssemblyLinePrediction] = []
    for key, indexed_modules in sorted(groups.items(), key=_group_sort_key):
        ordered, ordering_basis, ordering_confidence = _ordered_modules(
            record,
            indexed_modules,
            multi_cds_group=key[0] == "multi-cds",
        )
        modules = [module for _index, module in ordered]
        module_predictions = [
            _module_prediction(module, module_index, domains_for_module(record, module))
            for module_index, (_source_index, module) in enumerate(ordered, start=1)
        ]
        chain = tuple(call for module in module_predictions for call in module.monomer_calls)
        mass = _estimate_core_mass(tuple(module_predictions), chain)
        predictions.append(
            AssemblyLinePrediction(
                record_id=record.record_id,
                region_number=_region_number(record, modules),
                strand=_strand(modules),
                modules=tuple(module_predictions),
                chain=chain,
                ordering_basis=ordering_basis,
                ordering_confidence=ordering_confidence,
                caveats=(
                    *_prediction_caveats(
                        modules,
                        module_predictions,
                        ordering_basis=ordering_basis,
                        ordering_confidence=ordering_confidence,
                    ),
                    *_mass_caveats(mass),
                ),
                mass=mass,
            )
        )
    return predictions
