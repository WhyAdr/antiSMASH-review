"""Tabular export for evidence-only assembly-line predictions."""

from __future__ import annotations

import csv
import io

from antismash_review.assemblyline import (
    AssemblyLinePrediction,
    ModulePrediction,
    MonomerCall,
    predict_assembly_lines,
)
from antismash_review.models import Record

ASSEMBLYLINE_COLUMNS = (
    "record_id",
    "region_number",
    "assembly_line",
    "module_index",
    "locus_tags",
    "module_type",
    "complete",
    "starter",
    "final",
    "iterative",
    "multi_cds",
    "domain_names",
    "substrate",
    "monomer",
    "call_source",
    "call_confidence",
    "release_domains",
    "linear_core_mass_da",
    "head_to_tail_cyclic_candidate_mass_da",
    "mass_coverage",
    "mass_topology_assumption",
    "unresolved_components",
)


def _optional(value: object | None) -> object:
    return "" if value is None else value


def _list_cell(values: tuple[object, ...]) -> str:
    return "; ".join(str(value) for value in values)


def _row(
    prediction: AssemblyLinePrediction,
    module: ModulePrediction,
    call: MonomerCall | None,
    assembly_line: int,
) -> tuple[object, ...]:
    notes = (*module.warnings, *(call.notes if call else ()))
    mass = prediction.mass
    if mass is not None:
        notes = (*notes, *(f"unresolved monomer: {item}" for item in mass.unresolved_monomers))
    return (
        prediction.record_id,
        _optional(prediction.region_number),
        assembly_line,
        module.index,
        _list_cell(module.locus_tags),
        _optional(module.module_type),
        str(module.complete).lower(),
        str(module.starter).lower(),
        str(module.final).lower(),
        str(module.iterative).lower(),
        str(module.multi_cds).lower(),
        _list_cell(module.domain_names),
        _optional(call.substrate if call else None),
        _optional(call.monomer if call else None),
        _optional(call.source if call else None),
        _optional(call.confidence if call else None),
        _list_cell(module.release_domains),
        _optional(mass.linear_core_mass_da if mass else None),
        _optional(mass.head_to_tail_cyclic_candidate_mass_da if mass else None),
        _optional(mass.coverage_fraction if mass else None),
        _optional(mass.topology_assumption if mass else None),
        _list_cell(notes),
    )


def render_assemblyline_tsv(records: list[Record]) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    writer.writerow(ASSEMBLYLINE_COLUMNS)
    for record in records:
        for assembly_line, prediction in enumerate(predict_assembly_lines(record), start=1):
            for module in prediction.modules:
                calls = module.monomer_calls or (None,)
                for call in calls:
                    writer.writerow(_row(prediction, module, call, assembly_line))
    return output.getvalue()
