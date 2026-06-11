# avian-organoid-transcriptomics-reanalysis
Bioinformatics reanalysis of publicly available avian intestinal enteroid single-cell transcriptomic data (GSE283090) derived from broiler and layer chickens (Sun et al., 2025, Roslin Institute)
# Avian Organoid Transcriptomics Reanalysis
### Toward Breed-Relevant Models for East African Kienyeji Poultry

**Author:** Joy Manyasi Kabaka | ORCID: 0000-0002-2253-3970  
**Affiliation:** Independent Researcher, Nairobi, Kenya  
**Contact:** jmanyasa@gmail.com

---

## Project Overview

This repository presents a bioinformatics reanalysis of publicly available 
White Leghorn avian organoid transcriptomic datasets. The goal is to expose 
reference-genome bias and demonstrate the functional gaps that arise when 
commercial breed data is used as a universal proxy for all chicken breeds.

This work directly motivates a proposed extension: a variant-aware 
transcriptomic framework designed specifically for East African indigenous 
kienyeji breeds  genetic resources shaped by centuries of adaptation to 
thermal stress and endemic pathogens that commercial breeds cannot tolerate.

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
