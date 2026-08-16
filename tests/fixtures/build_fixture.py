from __future__ import annotations

import warnings
from pathlib import Path

from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, ExactPosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def _feature(
    feature_type: str,
    start: int,
    end: int,
    qualifiers: dict[str, list[str]] | None = None,
    *,
    fuzzy_start: bool = False,
) -> SeqFeature:
    left = BeforePosition(start) if fuzzy_start else ExactPosition(start)
    location = FeatureLocation(left, ExactPosition(end), strand=1)
    return SeqFeature(location=location, type=feature_type, qualifiers=qualifiers or {})


def write_synthetic_genbank(path: Path) -> Path:
    """Write a source-verified, biologically synthetic antiSMASH-like record."""
    record = SeqRecord(
        Seq("A" * 400),
        id="SYNTH.1",
        name="SYNTH.1",
        description="synthetic antiSMASH review fixture",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["structured_comment"] = {"antiSMASH-Data": {"Version": "8.0.4"}}
    record.features = [
        _feature("source", 0, 400, {"organism": ["synthetic bacterium"]}),
        _feature(
            "region",
            0,
            400,
            {
                "region_number": ["1"],
                "candidate_cluster_numbers": ["1", "1"],
                "product": ["NRPS", "synthetic product"],
                "contig_edge": ["true"],
            },
        ),
        _feature(
            "cand_cluster",
            20,
            380,
            {
                "candidate_cluster_number": ["1"],
                "protoclusters": ["1"],
                "product": ["NRPS"],
            },
        ),
        _feature(
            "protocluster",
            20,
            380,
            {
                "protocluster_number": ["1"],
                "candidate_cluster_numbers": ["1"],
            },
        ),
        _feature("proto_core", 0, 120, {"protocluster_number": ["1"]}),
        _feature(
            "CDS",
            0,
            90,
            {
                "locus_tag": ["SYN_CDS_1"],
                "gene_kind": ["biosynthetic"],
                "translation": ["M" * 30],
            },
            fuzzy_start=True,
        ),
        _feature(
            "gene",
            100,
            160,
            {
                "locus_tag": ["SYN_PSEUDO_1"],
                "pseudo": [""],
                "note": ["frameshift introduced in synthetic fixture"],
            },
        ),
        _feature(
            "aSDomain",
            30,
            80,
            {
                "domain_id": ["D1"],
                "aSDomain": ["A"],
                "aSTool": ["nrps_pks_domains"],
                "locus_tag": ["SYN_CDS_1"],
                "specificity": ["KR activity: inactive", "KR stereochemistry: C2"],
            },
        ),
        _feature(
            "aSDomain",
            80,
            120,
            {
                "domain_id": ["D1"],
                "aSDomain": ["TIGRFAM-domain"],
                "aSTool": ["tigrfam"],
                "locus_tag": ["SYN_CDS_1"],
            },
        ),
        _feature(
            "aSModule",
            30,
            120,
            {
                "domains": ["D1", "MISSING"],
                "locus_tags": ["SYN_CDS_1", "MISSING_CDS"],
                "complete": [""],
            },
        ),
        _feature(
            "CDS_motif",
            130,
            150,
            {
                "locus_tag": ["SYN_CDS_1"],
                "core_sequence": ["SYNCORE"],
            },
        ),
        _feature(
            "PFAM_domain",
            30,
            80,
            {
                "db_xref": ["PF00001.1"],
                "locus_tag": ["SYN_CDS_1"],
                "protein_start": ["1"],
                "protein_end": ["10"],
            },
        ),
        _feature(
            "PFAM_domain",
            30,
            80,
            {
                "db_xref": ["PF00001.2"],
                "locus_tag": ["SYN_CDS_1"],
                "protein_start": ["1"],
                "protein_end": ["10"],
            },
        ),
        _feature("tRNA", 170, 190),
        _feature("rRNA", 200, 220),
        _feature("repeat_region", 230, 250),
    ]
    # The synthetic fixture intentionally uses the long antiSMASH qualifier names;
    # Biopython warns while serializing them, but the parser must still preserve them.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Feature qualifier key '(candidate_cluster_numbers|candidate_cluster_number)' "
            r"is longer than maximum length specified by standard \(20 characters\)\.",
            category=BiopythonWarning,
        )
        SeqIO.write(record, path, "genbank")
    return path


def write_synthetic_cross_cds_monomer_genbank(path: Path) -> Path:
    """Write a synthetic GenBank record with a multi-CDS module containing duplicate pairings."""
    record = SeqRecord(
        Seq("A" * 600),
        id="SYNTH_CROSS.1",
        name="SYNTH_CROSS.1",
        description="synthetic cross-CDS duplicate monomer pairing fixture",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["structured_comment"] = {"antiSMASH-Data": {"Version": "8.0.4"}}
    record.features = [
        _feature("source", 0, 600, {"organism": ["synthetic bacterium"]}),
        _feature(
            "region",
            0,
            600,
            {
                "region_number": ["1"],
                "candidate_cluster_numbers": ["1"],
                "product": ["NRPS"],
            },
        ),
        _feature(
            "cand_cluster",
            0,
            600,
            {
                "candidate_cluster_number": ["1"],
                "protoclusters": ["1"],
                "product": ["NRPS"],
            },
        ),
        _feature(
            "protocluster",
            0,
            600,
            {
                "protocluster_number": ["1"],
                "candidate_cluster_numbers": ["1"],
            },
        ),
        _feature("proto_core", 0, 500, {"protocluster_number": ["1"]}),
        _feature(
            "CDS",
            50,
            250,
            {
                "locus_tag": ["CDS_A"],
                "gene_kind": ["biosynthetic"],
                "translation": ["M" * 66],
            },
        ),
        _feature(
            "CDS",
            300,
            500,
            {
                "locus_tag": ["CDS_B"],
                "gene_kind": ["biosynthetic"],
                "translation": ["M" * 66],
            },
        ),
        _feature(
            "aSDomain",
            60,
            200,
            {
                "domain_id": ["nrpspksdomains_CDS_A_A.1"],
                "aSDomain": ["AMP-binding"],
                "aSTool": ["nrps_pks_domains"],
                "locus_tag": ["CDS_A"],
                "specificity": ["substrate consensus: Orn"],
            },
        ),
        _feature(
            "aSDomain",
            320,
            450,
            {
                "domain_id": ["nrpspksdomains_CDS_B_PCP.1"],
                "aSDomain": ["PCP"],
                "aSTool": ["nrps_pks_domains"],
                "locus_tag": ["CDS_B"],
            },
        ),
        _feature(
            "aSModule",
            60,
            450,
            {
                "domains": [
                    "nrpspksdomains_CDS_A_A.1",
                    "nrpspksdomains_CDS_B_PCP.1",
                ],
                "locus_tags": ["CDS_A", "CDS_B"],
                "type": ["nrps"],
                "complete": [""],
                "monomer_pairings": ["Orn -> D-Orn", "Orn -> D-Orn"],
            },
        ),
    ]
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=BiopythonWarning,
        )
        SeqIO.write(record, path, "genbank")
    return path
