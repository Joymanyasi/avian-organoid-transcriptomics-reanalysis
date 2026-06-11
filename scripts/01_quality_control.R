# ============================================================
# SCRIPT 01: Quality Control
# Project: Avian Organoid Transcriptomics Reanalysis
# Dataset: GSE283090 (Sun et al., 2025, Roslin Institute)
# Author: Joy Manyasi Kabaka | ORCID: 0000-0002-2253-3970
# ============================================================

# --- INSTALL PACKAGES (run once) ---
# install.packages("BiocManager")
# BiocManager::install(c("Seurat", "ggplot2", "dplyr", "patchwork"))

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# --- LOAD DATA ---
# After downloading GSE283090_RAW.tar, extract and set your path here
# The tar file contains MTX and TSV files for each sample

# Example for one sample — repeat for all 8 samples
data_dir <- "data/GSM8655812/"  # Change to your extracted folder path

counts <- Read10X(data.dir = data_dir)

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "GSE283090",
  min.cells = 3,
  min.features = 200
)

# --- LABEL SAMPLE METADATA ---
seurat_obj$breed <- "Broiler"       # Change per sample: Broiler or Layer
seurat_obj$timepoint <- "D0"        # Change per sample: D0 or D3
seurat_obj$sample_id <- "GSM8655812" # Change per sample

# --- MITOCHONDRIAL QC ---
# High mitochondrial % = dying/poor quality cells
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj, 
  pattern = "^MT-"
)

# --- VISUALIZE QC METRICS ---
qc_plot <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

ggsave("results/01_QC_violin_plot.png", qc_plot, width = 12, height = 5)

# --- FILTER LOW QUALITY CELLS ---
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           percent.mt < 20
)

cat("Cells after QC filtering:", ncol(seurat_obj), "\n")

# --- SAVE FILTERED OBJECT ---
saveRDS(seurat_obj, "data/seurat_filtered.rds")

cat("Script 01 complete. Filtered object saved.\n")
