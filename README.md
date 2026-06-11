# avian-organoid-transcriptomics-reanalysis
Bioinformatics reanalysis of publicly available avian intestinal enteroid single-cell transcriptomic data (GSE283090) derived from broiler and layer chickens (Sun et al., 2025, Roslin Institute)
# Avian Organoid Transcriptomics Reanalysis
### Toward Breed-Relevant Models for East African Kienyeji Poultry

**Author:** Joy Manyasi Kabaka | ORCID: 0000-0002-2253-3970  
**Affiliation:** Independent Researcher, Nairobi, Kenya  
**Contact:** jmanyasa@gmail.com

---

## Project Overview
## Project Overview

The standard chicken reference genome (GRCg6a/GRCg7b) was derived from a 
single White Leghorn individual , a commercial inbred breed whose constrained 
genetic architecture does not represent the diversity of global chicken 
populations. Every avian transcriptomic study aligned against this reference 
inherits its bias, including the most comprehensive avian intestinal organoid 
atlas currently available (GSE283090, Sun et al., 2025, Roslin Institute).

This repository reanalyzes GSE283090 ,a single-cell RNA-seq atlas of 43,587 
broiler and layer intestinal enteroid cells through a variant-aware pipeline 
designed to expose what White Leghorn reference-genome bias conceals. By 
interrogating the ANP32A host-restriction locus and interferon-stimulated gene 
profiles, this reanalysis documents the functional gaps that arise when 
commercial breed data is treated as a universal proxy for all chickens.

This work directly motivates a proposed extension, a breed-specific 
transcriptomic framework for East African indigenous kienyeji chickens , 
genetic resources shaped by centuries of adaptation to thermal stress and 
endemic pathogens that commercial breeds cannot tolerate  establishing African 
data sovereignty as a foundational design principle.

---

## The Two-Chapter Story

### Chapter 1 — What Exists (This Repository)
- Reanalysis of public broiler and layer chicken intestinal enteroid scRNA-seq data (GSE283090)
- Standard alignment vs variant-aware pipeline comparison
- Interrogation of the ANP32A host-restriction locus
- Documenting where White Leghorn data falls short

### Chapter 2 — What Is Proposed (Future Work)
- Primary cell harvesting from kienyeji breeds
- iPSC reprogramming and 3D organoid derivation
- Population-specific SNP incorporation
- DSI governance aligned with CBD Nagoya Protocol

---

## Biological Question

Can existing public avian organoid transcriptomic data, reanalyzed 
through a variant-aware pipeline, reveal the functional limitations 
of White Leghorn reference genomes  and make the case for 
breed-specific African poultry models?

---

## Data Sources

| Dataset | Accession | Platform | Breed |
|---------|-----------|----------|-------|
| Sun et al. 2025 - Avian Intestinal Organoid scRNA-Seq Atlas | GSE283090 | scRNA-seq 10X Genomics | Broiler (Roslin) & Layer (Hy-Line Brown) |

- [x] Public dataset confirmed — GSE283090 (GEO, public April 2025)

---

## Repository Structure

data/           — Accession numbers and metadata (no raw files)
scripts/        — Analysis pipelines and code
results/        — Outputs, figures, tables
docs/           — Abstract, notes, governance documents

---

## Tools and Pipeline (Planned)

- FastQC — quality control
- STAR / HISAT2 — read alignment
- DESeq2 — differential expression
- GATK — variant-aware processing
- R / Python — downstream analysis

---

## DSI Governance Statement

This project embeds Digital Sequence Information (DSI) governance 
principles from inception. Any genomic assets derived from 
community-sourced kienyeji birds in future phases will remain 
legally anchored within African research ecosystems, consistent 
with the CBD Nagoya Protocol and principles of African data sovereignty.

---

## Status

- Conceptual framework established
- Public dataset — to be identified and confirmed
- Pipeline scripts — to be uploaded
- Reanalysis results — pending

---

## Citation

Kabaka, J.M. (2026). Beyond the White Leghorn: A Wet-to-Dry Organoid 
Framework for Functional Characterization of Climate-Resilient East 
African Kienyeji Poultry. Independent Research, Nairobi, Kenya.

---

## License
MIT License — open science, open access.
