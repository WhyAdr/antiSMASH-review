from __future__ import annotations

from .models import Diagnostic, Record, Severity


def _extend_consistency_diagnostics(
    record: Record,
    diagnostics: list[Diagnostic],
) -> None:
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


def review_record(record: Record) -> list[Diagnostic]:
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
