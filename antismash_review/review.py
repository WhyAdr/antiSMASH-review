from __future__ import annotations

from .models import Diagnostic, Record, Severity


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

    return diagnostics
