#!/usr/bin/env python3
"""
antiSMASH GenBank Parser & Secondary Metabolite Analyzer

Parses antiSMASH GenBank (.gbk) files using a hybrid engine:
  - Pre-parsing tag-aware line reflowing pass (fixes column ~79 continuation wrapping)
  - BioPython SeqIO integration (when installed) for robust feature & location parsing
  - Built-in zero-dependency fallback parser (when BioPython is absent)

Extracts 5 semantic layers of BGC data:
  1. Cluster Topology & Contig Edge Status
  2. Gene-level Inventory (Core vs Cargo / gene_kind tallies)
  3. NRPS/PKS Domain Architecture & Substrate Predictions
  4. Catalytic Domain Integrity (CDS_motif e-values & QC)
  5. Pfam Domain Annotations & Regulatory Synteny (LysR/HexS)
"""

import sys
import os
import re
import io
import json
import csv
import argparse
from pathlib import Path
from collections import defaultdict, Counter

# Try importing BioPython
HAS_BIOPYTHON = False
try:
    import Bio
    from Bio import SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

NO_SPACE_TAGS = {"translation", "domain_id", "locus_tag", "db_xref", "protein_id"}


def reflow_gbk_lines(lines):
    """
    Reflows GenBank flat-file continuation lines using tag-aware space re-insertion.
    Fixes column ~79 line wrapping without truncating tokens or inserting rogue spaces.
    """
    reflowed = []
    buf = ""
    curtag = ""

    for line in lines:
        line_str = line.rstrip("\r\n")

        # Feature start (5 spaces + feature key)
        if re.match(r"^     [A-Za-z_0-9]+", line_str):
            if buf:
                reflowed.append(buf)
            buf = line_str
            curtag = ""
        # Qualifier start (21 spaces + /tag= or /tag)
        elif re.match(r"^                     \/", line_str):
            if buf:
                reflowed.append(buf)
            buf = line_str
            m = re.match(r"^                     \/([A-Za-z_0-9]+)=", line_str)
            curtag = m.group(1) if m else ""
        # Continuation line (21 spaces + non-space/non-slash)
        elif re.match(r"^                     [^ \/]", line_str):
            content = line_str.lstrip()
            if curtag in NO_SPACE_TAGS:
                buf += content
            else:
                buf += " " + content
        else:
            if buf:
                reflowed.append(buf)
                buf = ""
            reflowed.append(line_str)

    if buf:
        reflowed.append(buf)
    return reflowed


def parse_qualifiers_raw(qual_lines):
    """Fallback parser for qualifier lines."""
    quals = defaultdict(list)
    for qstr in qual_lines:
        qstr_clean = qstr.strip()
        if not qstr_clean.startswith("/"):
            continue
        if "=" in qstr_clean:
            key, val = qstr_clean[1:].split("=", 1)
            val = val.strip('"')
            quals[key].append(val)
        else:
            key = qstr_clean[1:]
            quals[key].append("True")
    return quals


def parse_gbk_features_fallback(reflowed_lines):
    """Fallback parser when BioPython is not installed."""
    features = []
    in_features = False
    current_feature = None

    for line in reflowed_lines:
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN") or line.startswith("//"):
            in_features = False
            if current_feature:
                features.append(current_feature)
                current_feature = None
            break

        if not in_features:
            continue

        m_feat = re.match(r"^     ([A-Za-z_0-9]+)\s+(.+)", line)
        if m_feat:
            if current_feature:
                features.append(current_feature)
            key = m_feat.group(1).strip()
            loc = m_feat.group(2).strip()
            current_feature = {
                "key": key,
                "location": loc,
                "qual_lines": []
            }
        elif line.startswith("                     /") and current_feature:
            current_feature["qual_lines"].append(line)

    if current_feature:
        features.append(current_feature)

    for feat in features:
        feat["qualifiers"] = parse_qualifiers_raw(feat.pop("qual_lines", []))

    return features


def parse_gbk_features_biopython(reflowed_lines):
    """Parses reflowed GenBank records using BioPython SeqIO."""
    reflowed_text = "\n".join(reflowed_lines)
    records = list(SeqIO.parse(io.StringIO(reflowed_text), "genbank"))
    
    features = []
    for record in records:
        for feat in record.features:
            # Convert BioPython qualifiers dict
            quals = {k: v if isinstance(v, list) else [v] for k, v in feat.qualifiers.items()}
            features.append({
                "key": feat.type,
                "location": str(feat.location),
                "qualifiers": quals
            })
    return features


def analyze_antismash_record(filepath):
    """Analyzes an antiSMASH GenBank file across all 5 semantic layers."""
    path = Path(filepath)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        raw_lines = f.readlines()

    # Pre-parsing reflowing pass
    reflowed = reflow_gbk_lines(raw_lines)

    # Use BioPython if available, otherwise fallback
    engine_used = "BioPython" if HAS_BIOPYTHON else "Built-in Fallback"
    if HAS_BIOPYTHON:
        try:
            features = parse_gbk_features_biopython(reflowed)
        except Exception:
            features = parse_gbk_features_fallback(reflowed)
            engine_used = "Built-in Fallback (BioPython parse exception)"
    else:
        features = parse_gbk_features_fallback(reflowed)

    analysis = {
        "filename": path.name,
        "filepath": str(path.resolve()),
        "parser_engine": engine_used,
        "cluster_topology": [],
        "gene_kinds": Counter(),
        "genes_and_cds": [],
        "nrps_pks_domains": [],
        "cds_motifs": [],
        "pfam_domains": [],
        "tfbs_features": [],
        "lysr_regulators": [],
        "qc_alerts": []
    }

    for feat in features:
        key = feat["key"]
        loc = feat["location"]
        quals = feat["qualifiers"]

        # 1. Cluster Topology
        if key in ("region", "cand_cluster", "protocluster", "proto_core"):
            prod = quals.get("product", ["Unknown"])[0]
            contig_edge = quals.get("contig_edge", ["Unknown"])[0]
            region_num = quals.get("region_number", ["1"])[0]
            rules = quals.get("rules", [])
            core_loc = quals.get("core_location", [loc])[0]

            topo_entry = {
                "layer": key,
                "product": prod,
                "contig_edge": contig_edge,
                "region_number": region_num,
                "location": loc,
                "core_location": core_loc,
                "rules": rules
            }
            analysis["cluster_topology"].append(topo_entry)

            if str(contig_edge).lower() == "true":
                analysis["qc_alerts"].append(
                    f"CRITICAL: Region {region_num} ({prod}) is flagged contig_edge='True' (truncated at assembly boundary)."
                )

        # 2. Gene-Level Cross-Walk (gene / CDS)
        elif key == "CDS":
            locus = quals.get("locus_tag", [""])[0]
            gene = quals.get("gene", [""])[0]
            product = quals.get("product", [""])[0]
            kind = quals.get("gene_kind", ["other"])[0]
            inference = quals.get("inference", [])

            analysis["gene_kinds"][kind] += 1

            cds_entry = {
                "locus_tag": locus,
                "gene": gene,
                "product": product,
                "gene_kind": kind,
                "location": loc,
                "inference": inference
            }
            analysis["genes_and_cds"].append(cds_entry)

            prod_lower = product.lower()
            if "lysr" in prod_lower or "auto-aggregation" in prod_lower or "hexs" in prod_lower or gene.lower() == "hexs":
                analysis["lysr_regulators"].append({
                    "locus_tag": locus,
                    "gene": gene,
                    "product": product,
                    "location": loc
                })

        # 3. NRPS/PKS Domain Architecture (aSDomain / aSModule)
        elif key in ("aSDomain", "aSModule"):
            domain_id = quals.get("domain_id", [""])[0]
            asdomain = quals.get("aSDomain", [""])[0]
            subtype = quals.get("domain_subtype", [""])[0]
            spec = quals.get("specificity", [])
            monomers = quals.get("monomer_pairings", [])
            smiles = quals.get("SMILES", [])
            locus = quals.get("locus_tag", [""])[0]

            dom_entry = {
                "layer": key,
                "domain_id": domain_id,
                "domain": asdomain or subtype,
                "specificity": spec,
                "monomer_pairings": monomers,
                "smiles": smiles,
                "locus_tag": locus,
                "location": loc
            }
            analysis["nrps_pks_domains"].append(dom_entry)

        # 4. Catalytic Motifs (CDS_motif)
        elif key == "CDS_motif":
            motif = quals.get("motif", quals.get("label", [""])[0])[0]
            evalue = quals.get("evalue", quals.get("evalue_label", [""])[0])[0]
            score = quals.get("score", [""])[0]
            locus = quals.get("locus_tag", [""])[0]

            motif_entry = {
                "motif": motif,
                "evalue": evalue,
                "score": score,
                "locus_tag": locus,
                "location": loc
            }
            analysis["cds_motifs"].append(motif_entry)

            try:
                if evalue and float(evalue) > 1e-3:
                    analysis["qc_alerts"].append(
                        f"WARNING: Low-confidence catalytic motif '{motif}' (e-value: {evalue}) in locus {locus}."
                    )
            except ValueError:
                pass

        # 5. Pfam Domains
        elif key == "PFAM_domain":
            db_xref = quals.get("db_xref", [])
            desc = quals.get("description", [""])[0]
            evalue = quals.get("evalue", [""])[0]
            locus = quals.get("locus_tag", [""])[0]

            analysis["pfam_domains"].append({
                "db_xref": db_xref,
                "description": desc,
                "evalue": evalue,
                "locus_tag": locus,
                "location": loc
            })

        # TFBS Hits
        elif key == "misc_feature":
            note = quals.get("note", [""])[0]
            if "TFBS" in note or "binding" in note.lower() or "tfbs" in note.lower():
                analysis["tfbs_features"].append({
                    "note": note,
                    "location": loc
                })

    return analysis


def generate_markdown_summary(analysis):
    """Generates a structured human-readable markdown summary report."""
    md = []
    md.append(f"# antiSMASH Analysis Report: `{analysis['filename']}`")
    md.append(f"*Parser Engine Used*: `{analysis['parser_engine']}`\n")

    if analysis["qc_alerts"]:
        md.append("## 🚨 Quality Control Alerts & Flags")
        for alert in analysis["qc_alerts"]:
            md.append(f"- {alert}")
        md.append("")

    md.append("## 1. Cluster Topology & Detection")
    if analysis["cluster_topology"]:
        for topo in analysis["cluster_topology"]:
            md.append(f"- **Layer**: `{topo['layer']}` | **Product**: `{topo['product']}` | **Contig Edge**: `{topo['contig_edge']}`")
            md.append(f"  - Location: `{topo['location']}` (Core: `{topo['core_location']}`)")
            if topo["rules"]:
                md.append(f"  - Detection Rules: {', '.join(topo['rules'])}")
    else:
        md.append("No explicit region/protocluster records found.")
    md.append("")

    md.append("## 2. Gene-Level Functional Classification (`gene_kind`)")
    md.append("| Gene Kind | Count | Description |")
    md.append("| --- | --- | --- |")
    kind_desc = {
        "biosynthetic": "Core biosynthetic gene(s)",
        "biosynthetic-additional": "Additional biosynthetic / tailoring gene(s)",
        "regulatory": "Pathway-specific regulatory gene(s)",
        "transport": "Transporter / efflux gene(s)",
        "other": "Other cargo / flanking gene(s)"
    }
    for kind, count in analysis["gene_kinds"].items():
        md.append(f"| `{kind}` | {count} | {kind_desc.get(kind, 'Unclassified')} |")
    md.append("")

    md.append("## 3. Regulatory Synteny & Candidate Regulators (LysR/HexS)")
    if analysis["lysr_regulators"]:
        for reg in analysis["lysr_regulators"]:
            md.append(f"- **Locus**: `{reg['locus_tag']}` | **Gene**: `{reg['gene'] or 'N/A'}` | **Product**: `{reg['product']}`")
            md.append(f"  - Location: `{reg['location']}`")
    else:
        md.append("No LysR-family / HexS ortholog candidates detected in immediate BGC vicinity.")
    md.append("")

    md.append("## 4. NRPS/PKS Domain Architecture & Specificity")
    if analysis["nrps_pks_domains"]:
        for dom in analysis["nrps_pks_domains"]:
            spec_str = "; ".join(dom["specificity"]) if dom["specificity"] else "N/A"
            mon_str = "; ".join(dom["monomer_pairings"]) if dom["monomer_pairings"] else "N/A"
            md.append(f"- **Domain**: `{dom['domain']}` | **Locus**: `{dom['locus_tag']}`")
            md.append(f"  - Specificity: `{spec_str}` | Monomers: `{mon_str}`")
    else:
        md.append("No modular NRPS/PKS domain annotations (`aSDomain`) found.")
    md.append("")

    md.append("## 5. Catalytic Domain Integrity (`CDS_motif`)")
    md.append(f"Total conserved catalytic motifs detected: **{len(analysis['cds_motifs'])}**")
    if analysis["cds_motifs"]:
        motif_names = set(m["motif"] for m in analysis["cds_motifs"])
        md.append(f"Motifs present: {', '.join(sorted(motif_names))}")
    md.append("")

    md.append("## 6. Pfam Domains & TFBS Notes")
    md.append(f"- Pfam domains annotated across region: **{len(analysis['pfam_domains'])}**")
    md.append(f"- Predicted TFBS hits (`misc_feature`): **{len(analysis['tfbs_features'])}**")
    if analysis["tfbs_features"]:
        md.append("  - *Note*: antiSMASH TFBS Finder uses LogoMotif (Actinobacteria-focused). Exercise caution when interpreting TFBS predictions in Enterobacterales/Serratia.")

    return "\n".join(md)


def main():
    parser = argparse.ArgumentParser(
        description="antiSMASH GenBank Flat-File Parser & Secondary Metabolite Analyzer (Hybrid Engine)"
    )
    parser.add_argument("input_path", help="Path to input .gbk file or directory containing .gbk files")
    parser.add_argument(
        "--format",
        choices=["summary", "json", "tsv"],
        default="summary",
        help="Output format (summary, json, tsv)"
    )
    parser.add_argument("--output", help="Output file path (default: stdout)")

    args = parser.parse_args()

    input_path = Path(args.input_path)
    if not input_path.exists():
        print(f"Error: Input path '{args.input_path}' does not exist.", file=sys.stderr)
        sys.exit(1)

    files_to_parse = []
    if input_path.is_file():
        files_to_parse.append(input_path)
    else:
        files_to_parse.extend(sorted(input_path.glob("*.gbk")))
        files_to_parse.extend(sorted(input_path.glob("*.gb")))

    if not files_to_parse:
        print(f"Error: No .gbk or .gb files found at '{args.input_path}'.", file=sys.stderr)
        sys.exit(1)

    results = [analyze_antismash_record(f) for f in files_to_parse]

    out_text = ""
    if args.format == "json":
        for r in results:
            r["gene_kinds"] = dict(r["gene_kinds"])
        out_text = json.dumps(results if len(results) > 1 else results[0], indent=2)

    elif args.format == "summary":
        summaries = [generate_markdown_summary(r) for r in results]
        out_text = "\n\n---\n\n".join(summaries)

    elif args.format == "tsv":
        lines = ["filename\tengine\tregion_product\tcontig_edge\tcore_genes\ttotal_genes\tnrps_pks_domains\tlysr_regulators\tqc_alerts"]
        for r in results:
            fname = r["filename"]
            engine = r["parser_engine"]
            prod = r["cluster_topology"][0]["product"] if r["cluster_topology"] else "Unknown"
            edge = r["cluster_topology"][0]["contig_edge"] if r["cluster_topology"] else "Unknown"
            core = r["gene_kinds"].get("biosynthetic", 0)
            total = sum(r["gene_kinds"].values())
            doms = len(r["nrps_pks_domains"])
            lysr = len(r["lysr_regulators"])
            alerts = "; ".join(r["qc_alerts"]) if r["qc_alerts"] else "None"
            lines.append(f"{fname}\t{engine}\t{prod}\t{edge}\t{core}\t{total}\t{doms}\t{lysr}\t{alerts}")
        out_text = "\n".join(lines)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(out_text)
        print(f"Analysis written to '{args.output}'.")
    else:
        print(out_text)


if __name__ == "__main__":
    main()
