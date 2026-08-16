"""Human-readable export for evidence-only assembly-line predictions."""

from __future__ import annotations

from antismash_review.assemblyline import predict_assembly_lines
from antismash_review.models import Record


def render_assemblyline_markdown(records: list[Record]) -> str:
    lines = ["# antiSMASH assembly-line evidence", ""]
    for record in records:
        lines.extend([f"## `{record.record_id}`", ""])
        predictions = predict_assembly_lines(record)
        if not predictions:
            lines.extend(["- No antiSMASH modules were parsed.", ""])
            continue
        for line_index, prediction in enumerate(predictions, start=1):
            region = (
                str(prediction.region_number)
                if prediction.region_number is not None
                else "unresolved"
            )
            lines.extend(
                [
                    f"### Assembly line {line_index} (region {region})",
                    "",
                    f"- Ordering: `{prediction.ordering_basis}` "
                    f"({prediction.ordering_confidence} confidence)",
                    "- Strand: "
                    f"`{prediction.strand if prediction.strand is not None else 'unresolved'}`",
                    (
                        f"- Mass coverage: `{prediction.mass.coverage_fraction:.3f}`; "
                        f"linear candidate: `{prediction.mass.linear_core_mass_da}`; "
                        f"cyclic candidate: `"
                        f"{prediction.mass.head_to_tail_cyclic_candidate_mass_da}`"
                        if prediction.mass is not None
                        else "- Mass: not calculated"
                    ),
                    (
                        "| Module | Locus tags | Type | Raw calls | Incorporation | "
                        "Pairing status | Integrity flags | Domains | Release evidence |"
                    ),
                    "|---:|---|---|---|---|---|---|---|---|",
                ]
            )
            for module in prediction.modules:
                raw_calls = ", ".join(call.display for call in module.monomer_calls) or "?"
                incorporation = module.incorporation_call.display or "?"
                pairing_status = module.pairing_status
                integrity = ", ".join(module.integrity_flags) or "none"
                domains = ", ".join(module.domain_names) or "none resolved"
                release = ", ".join(module.release_domains) or "none"
                lines.append(
                    f"| {module.index} | {', '.join(module.locus_tags) or 'unassigned'} | "
                    f"{module.module_type or 'unknown'} | {raw_calls} | {incorporation} | "
                    f"{pairing_status} | {integrity} | {domains} | {release} |"
                )
            if prediction.caveats:
                lines.extend(["", "Caveats:"])
                lines.extend(f"- {caveat}" for caveat in prediction.caveats)
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"
