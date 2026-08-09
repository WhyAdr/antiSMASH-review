"""Supported public API for parsing, reviewing, and exporting antiSMASH GenBank records."""

from ._version import __version__
from .architecture import assess_architecture
from .assemblyline import predict_assembly_lines
from .cohort import build_cohort
from .exporters.bed import render_bed
from .exporters.entity_tables import (
    render_clusterblast_tsv,
    render_domain_tsv,
    render_gene_tsv,
)
from .exporters.gff3 import render_gff3
from .exporters.json_export import dumps_records
from .exporters.markdown import render_records
from .exporters.provenance import dumps_provenance as render_provenance_json
from .exporters.tables import render_tsv
from .genbank import GenBankParseError, parse_genbank
from .models import Diagnostic, Record, Severity
from .review import ReviewFinding, review_findings, review_record

__all__ = [
    "Diagnostic",
    "GenBankParseError",
    "Record",
    "ReviewFinding",
    "Severity",
    "__version__",
    "assess_architecture",
    "build_cohort",
    "dumps_records",
    "parse_genbank",
    "predict_assembly_lines",
    "render_clusterblast_tsv",
    "render_bed",
    "render_domain_tsv",
    "render_gene_tsv",
    "render_gff3",
    "render_provenance_json",
    "render_records",
    "render_tsv",
    "review_findings",
    "review_record",
]
