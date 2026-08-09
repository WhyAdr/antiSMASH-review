from __future__ import annotations

from dataclasses import dataclass

from .architecture import assess_architecture
from .locations import modules_in_collection, overlaps
from .models import Diagnostic, Location, RawFeature, Record, Severity


@dataclass(slots=True, frozen=True)
class ReviewFinding:
    """A diagnostic plus structured entity context for browser-track export."""

    diagnostic: Diagnostic
    location: Location | None
    entity_type: str | None
    entity_id: str | None


def _locations_overlap(left: Location, right: Location) -> bool:
    return any(
        left_part.start < right_part.end and right_part.start < left_part.end
        for left_part in left.parts
        for right_part in right.parts
    )


def _raw_label(raw: RawFeature) -> str:
    for key in ("locus_tag", "gene"):
        values = raw.qualifiers.get(key, ())
        if values and values[0]:
            return values[0]
    return "unlabelled gene"


def _extend_consistency_diagnostics(
    record: Record,
    diagnostics: list[Diagnostic],
) -> None:
    for raw in record.raw_features:
        if raw.feature_type != "gene" or "pseudo" not in raw.qualifiers:
            continue
        overlapping_regions = [
            region for region in record.regions if _locations_overlap(raw.location, region.location)
        ]
        if not overlapping_regions:
            continue
        region_numbers = sorted(
            region.number for region in overlapping_regions if region.number is not None
        )
        region_text = ", ".join(str(number) for number in region_numbers) or "unnumbered"
        diagnostics.append(
            Diagnostic(
                code="pseudogene_in_cluster",
                severity=Severity.WARNING,
                message=(
                    f"Pseudo gene {_raw_label(raw)} overlaps antiSMASH region(s): "
                    f"{region_text}; inspect frameshift or annotation evidence"
                ),
                source=str(record.source_path),
                record_id=record.record_id,
                feature_index=raw.feature_index,
            )
        )

    gene_locus_tags = {gene.locus_tag for gene in record.genes if gene.locus_tag}
    for module_index, module in enumerate(record.modules):
        orphan_tags = sorted({tag for tag in module.locus_tags if tag not in gene_locus_tags})
        if orphan_tags:
            diagnostics.append(
                Diagnostic(
                    code="orphan_module_locus",
                    severity=Severity.WARNING,
                    message=(
                        f"Module {module_index + 1} references locus tags absent from "
                        f"the CDS set: {', '.join(orphan_tags)}"
                    ),
                    source=str(record.source_path),
                    record_id=record.record_id,
                )
            )

    architecture_products = sorted(
        {
            product
            for region in record.regions
            for product in region.products
            if "nrps" in product.casefold() or "pks" in product.casefold()
        }
    )
    if architecture_products and not record.nrps_pks_domains:
        diagnostics.append(
            Diagnostic(
                code="missing_nrps_pks_architecture",
                severity=Severity.WARNING,
                message=(
                    "Region products imply NRPS/PKS architecture "
                    f"({', '.join(architecture_products)}), but no domains from "
                    "aSTool=nrps_pks_domains were parsed"
                ),
                source=str(record.source_path),
                record_id=record.record_id,
            )
        )

    for assessment in assess_architecture(record):
        if not assessment.missing_slots:
            continue
        observed = ", ".join(assessment.observed_slots) or "none"
        missing = ", ".join(assessment.missing_slots)
        diagnostics.append(
            Diagnostic(
                code="architecture_core_domain_missing",
                severity=Severity.WARNING,
                message=(
                    f"{assessment.product} {assessment.scope} lacks expected core domain "
                    f"slot(s): {missing}. Observed slots: {observed}. This is an "
                    "architecture-consistency warning, not evidence that the locus is "
                    "nonfunctional."
                ),
                source=str(record.source_path),
                record_id=record.record_id,
            )
        )


def _review_record_diagnostics(record: Record) -> list[Diagnostic]:
    """Apply conservative, evidence-scoped review rules to one record."""
    diagnostics = list(record.diagnostics)

    edge_regions = [region for region in record.regions if region.contig_edge is True]
    if edge_regions:
        diagnostics.append(
            Diagnostic(
                code="context_reaches_record_edge",
                severity=Severity.NOTICE,
                message="BGC context reaches a record boundary; gene context may be incomplete",
                source=str(record.source_path),
                record_id=record.record_id,
            )
        )

    for core in record.proto_cores:
        if core.location.start == 0 or core.location.end == record.length:
            diagnostics.append(
                Diagnostic(
                    code="core_reaches_record_edge",
                    severity=Severity.WARNING,
                    message="A proto-core reaches a record boundary",
                    source=str(record.source_path),
                    record_id=record.record_id,
                )
            )

    for gene in record.genes:
        at_edge = gene.location.start == 0 or gene.location.end == record.length
        if at_edge and gene.location.partial:
            diagnostics.append(
                Diagnostic(
                    code="partial_cds_at_edge",
                    severity=Severity.WARNING,
                    message=f"Partial CDS reaches a boundary: {gene.locus_tag or 'unlabelled CDS'}",
                    source=str(record.source_path),
                    record_id=record.record_id,
                )
            )

    _extend_consistency_diagnostics(record, diagnostics)

    return diagnostics


def _assessment_entity(
    record: Record,
    assessment: object,
) -> tuple[Location | None, str | None, str | None]:
    region_number = getattr(assessment, "region_number", None)
    module_index = getattr(assessment, "module_index", None)
    regions = [
        region
        for region in record.regions
        if region_number is None or region.number == region_number
    ]
    if not regions:
        return None, None, None
    region = regions[0]
    if module_index is not None:
        modules = modules_in_collection(record, region)
        if 0 < module_index <= len(modules):
            module = modules[module_index - 1]
            return module.location, "aSModule", module.locus_tags[0] if module.locus_tags else None
    return region.location, "region", str(region.number) if region.number is not None else None


def _finding_candidates(
    record: Record,
) -> dict[str, list[tuple[Location | None, str | None, str | None]]]:
    """Build structural locations in the same append order as review diagnostics."""

    candidates: dict[str, list[tuple[Location | None, str | None, str | None]]] = {}

    def add(
        code: str,
        location: Location | None,
        entity_type: str | None,
        entity_id: str | None,
    ) -> None:
        candidates.setdefault(code, []).append((location, entity_type, entity_id))

    edge_regions = [region for region in record.regions if region.contig_edge is True]
    if edge_regions:
        region = edge_regions[0]
        add("context_reaches_record_edge", region.location, "region", str(region.number))
    for core in record.proto_cores:
        if core.location.start == 0 or core.location.end == record.length:
            add("core_reaches_record_edge", core.location, "proto_core", str(core.number))
    for gene in record.genes:
        at_edge = gene.location.start == 0 or gene.location.end == record.length
        if at_edge and gene.location.partial:
            add("partial_cds_at_edge", gene.location, "CDS", gene.locus_tag)

    for raw in record.raw_features:
        if raw.feature_type != "gene" or "pseudo" not in raw.qualifiers:
            continue
        if any(overlaps(raw.location, region.location) for region in record.regions):
            add("pseudogene_in_cluster", raw.location, "gene", _raw_label(raw))

    gene_locus_tags = {gene.locus_tag for gene in record.genes if gene.locus_tag}
    for module in record.modules:
        orphan_tags = {tag for tag in module.locus_tags if tag not in gene_locus_tags}
        if orphan_tags:
            add(
                "orphan_module_locus",
                module.location,
                "aSModule",
                module.locus_tags[0] if module.locus_tags else None,
            )

    architecture_products = {
        product
        for region in record.regions
        for product in region.products
        if "nrps" in product.casefold() or "pks" in product.casefold()
    }
    if architecture_products and not record.nrps_pks_domains:
        matched_region = next(
            (
                region
                for region in record.regions
                if any(product in architecture_products for product in region.products)
            ),
            None,
        )
        add(
            "missing_nrps_pks_architecture",
            matched_region.location if matched_region is not None else None,
            "region" if matched_region is not None else None,
            str(matched_region.number)
            if matched_region is not None and matched_region.number is not None
            else None,
        )

    for assessment in assess_architecture(record):
        if assessment.missing_slots:
            add("architecture_core_domain_missing", *_assessment_entity(record, assessment))
    return candidates


def review_findings(record: Record) -> list[ReviewFinding]:
    """Return review diagnostics with entity locations where they are known."""

    diagnostics = _review_record_diagnostics(record)
    candidates = _finding_candidates(record)
    findings: list[ReviewFinding] = []
    source_diagnostic_count = len(record.diagnostics)
    for index, diagnostic in enumerate(diagnostics):
        if index < source_diagnostic_count:
            location, entity_type, entity_id = None, None, None
        else:
            options = candidates.get(diagnostic.code, [])
            location, entity_type, entity_id = options.pop(0) if options else (None, None, None)
        findings.append(
            ReviewFinding(
                diagnostic=diagnostic,
                location=location,
                entity_type=entity_type,
                entity_id=entity_id,
            )
        )
    return findings


def review_record(record: Record) -> list[Diagnostic]:
    """Apply conservative, evidence-scoped review rules to one record."""

    return [finding.diagnostic for finding in review_findings(record)]
