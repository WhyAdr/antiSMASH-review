---
name: antismash-genbank-parser
description: Specialized parser and analysis tool for antiSMASH GenBank (.gbk) files. Parses cluster topology, contig_edge status, gene_kind classifications, NRPS/PKS domain architecture, catalytic CDS_motifs, and Pfam domains into structured summaries.
---

# antismash-genbank-parser

A specialized parser and analytical workflow for antiSMASH GenBank (.gbk) flat-files. Designed for secondary metabolite workflows, BGC panel audits, and biosynthetic gene cluster dissection.

## When to Use

- When parsing and analyzing antiSMASH GenBank output files (`.gbk` / `.gb`).
- When performing batch audits of Biosynthetic Gene Clusters (BGCs) across multiple isolate genomes.
- When inspecting NRPS/PKS domain architecture, monomer predictions, and catalytic domain integrity.
- When tallying core vs. cargo genes using antiSMASH `gene_kind` classifications.
- When evaluating BGC truncation status (`contig_edge="True"`) or searching for adjacent regulatory genes (e.g. LysR-family / `hexS` orthologs).

## Hybrid Engine Architecture

The bundled script `scripts/parse_antismash_gbk.py` uses a **hybrid parsing engine**:
1. **Tag-Aware Line Reflowing (Pre-pass)**: Fixes column ~79 continuation line wrapping using specific re-insertion rules before passing to feature parsing.
2. **BioPython Integration (Preferred)**: Uses `Bio.SeqIO` to parse feature objects, coordinates, strands, and qualifiers if BioPython is installed.
3. **Built-in Fallback**: Falls back gracefully to a standard-library parser if BioPython is unavailable in the environment.

### Core Line Reflowing Rules
1. **Free-text qualifiers** (`/product=`, `/note=`, `/NRPS_PKS=`): Wrap at space boundaries. Re-insert a space when rejoining continuation lines.
2. **Fixed token qualifiers** (`/translation=`, `/domain_id=`, `/locus_tag=`, `/db_xref=`): Wrap mid-token. Rejoin directly **without** adding spaces.

## 5 Semantic Layers Parsed

1. **Cluster Topology (`region`, `cand_cluster`, `protocluster`, `proto_core`)**
   - Extracts product type, detection rules, core location, and `contig_edge` flag.
   - **BGC Completeness Check**: A BGC with `contig_edge="True"` is truncated at an assembly boundary and must be cross-checked against contig breakpoints before assuming gene content is complete.

2. **Gene-Level Cross-Walk (`gene`, `CDS`)**
   - Cross-references locus tags, gene symbols, products, EC numbers, and Bakta evidence codes.
   - Groups CDS features by antiSMASH `gene_kind` calls (`biosynthetic`, `biosynthetic-additional`, `regulatory`, `transport`, `other`) to provide instant "core vs. cargo" counts.

3. **NRPS/PKS Domain Architecture (`aSDomain`, `aSModule`)**
   - Extracts domain subtypes (e.g., `Condensation_Starter`, `AMP-binding`, `PCP`, `Thioesterase`), Stachelhaus code specificity predictions, `monomer_pairings` (e.g., `Ser -> Ser`), and SMILES structures.
   - **Sanity Check**: Verifies predicted substrate specificity against characterized chemical structures (e.g., confirming `Ser` calls for `swrW` serrawettin W1 synthetase).

4. **Catalytic Domain QC (`CDS_motif`)**
   - Extracts 11 conserved catalytic motifs (C2/C3/C5/C67, NRPS-A_a2/a3/a6/a8, NRPS-E2) and e-values.
   - **Integrity Alert**: Missing or high-e-value motifs (>1e-3) inside an otherwise "complete" module indicate potential frameshifts or active-site degeneracy.

5. **Broad Annotations & Regulatory Context (`PFAM_domain`, `misc_feature`)**
   - **Pfam Domains**: Captures Pfam-A hits across the full extraction window.
   - **TFBS Caveat**: antiSMASH TFBS Finder uses LogoMotif (Actinobacteria-focused). TFBS predictions in Enterobacterales (*Serratia*, *E. coli*) should be treated as hypothesis-generating.
   - **Regulatory Synteny**: Scans BGC flanks for local transcriptional regulators (e.g. LysR-family / `hexS` orthologs).

## CLI Script Usage

```bash
# Generate human-readable Markdown report for a single BGC region
python3 scripts/parse_antismash_gbk.py region011.gbk --format summary

# Batch parse an entire directory of BGC regions into a TSV table
python3 scripts/parse_antismash_gbk.py /path/to/antiSMASH_results/ --format tsv --output bgc_summary_table.tsv

# Export full structured data as JSON
python3 scripts/parse_antismash_gbk.py region011.gbk --format json --output region011.json
```
