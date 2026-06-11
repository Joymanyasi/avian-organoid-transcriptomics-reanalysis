# ============================================================
# SCRIPT 01: Load All 8 Samples, Merge & Quality Control
# Project: Avian Organoid Transcriptomics Reanalysis
# Dataset: GSE283090 (Sun et al., 2025, Roslin Institute)
# Author: Joy Manyasi Kabaka | ORCID: 0000-0002-2253-3970
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# --- SAMPLE METADATA ---
# Update data_root to wherever you extracted GSE283090_RAW.tar
data_root <- "data/raw/"

sample_info <- data.frame(
  sample_id  = c("GSM8655812","GSM8655813","GSM8655814","GSM8655815",
                 "GSM8655816","GSM8655817","GSM8655818","GSM8655819"),
  breed      = c("Broiler","Broiler","Layer","Layer",
                 "Broiler","Broiler","Layer","Layer"),
  timepoint  = c("D0","D0","D0","D0","D3","D3","D3","D3"),
  replicate  = c("Rep1","Rep2","Rep1","Rep2","Rep1","Rep2","Rep1","Rep2"),
  stringsAsFactors = FALSE
)

# --- LOAD AND CREATE SEURAT OBJECTS ---
seurat_list <- list()

for (i in seq_len(nrow(sample_info))) {
  sid  <- sample_info$sample_id[i]
  path <- file.path(data_root, sid)
  
  cat("Loading sample:", sid, "\n")
  
  counts <- Read10X(data.dir = path)
  
  obj <- CreateSeuratObject(
    counts       = counts,
    project      = sid,
    min.cells    = 3,
    min.features = 200
  )
  
  obj$sample_id  <- sid
  obj$breed      <- sample_info$breed[i]
  obj$timepoint  <- sample_info$timepoint[i]
  obj$replicate  <- sample_info$replicate[i]
  
  seurat_list[[sid]] <- obj
}

# --- MERGE ALL 8 SAMPLES INTO ONE OBJECT ---
merged <- merge(
  x            = seurat_list[[1]],
  y            = seurat_list[-1],
  add.cell.ids = sample_info$sample_id,
  project      = "GSE283090"
)

cat("Total cells after merge:", ncol(merged), "\n")

# --- MITOCHONDRIAL CONTENT ---
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

# --- VISUALIZE QC BEFORE FILTERING ---
p_before <- VlnPlot(
  merged,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol     = 3,
  pt.size  = 0
) + plot_annotation(title = "QC Metrics Before Filtering — All 8 Samples")

ggsave("results/01_QC_before_filtering.png", p_before, width = 16, height = 5)

# --- FILTER LOW QUALITY CELLS ---
merged_filtered <- subset(
  merged,
  subset = nFeature_RNA > 200  &
           nFeature_RNA < 6000 &
           nCount_RNA   > 500  &
           percent.mt   < 20
)

cat("Cells before filtering:", ncol(merged), "\n")
cat("Cells after  filtering:", ncol(merged_filtered), "\n")
cat("Cells removed:         ", ncol(merged) - ncol(merged_filtered), "\n")

# --- VISUALIZE QC AFTER FILTERING ---
p_after <- VlnPlot(
  merged_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol     = 3,
  pt.size  = 0
) + plot_annotation(title = "QC Metrics After Filtering — All 8 Samples")

ggsave("results/01_QC_after_filtering.png", p_after, width = 16, height = 5)

# --- SAVE ---
saveRDS(merged_filtered, "data/seurat_merged_filtered.rds")
cat("Script 01 complete. Merged filtered object saved.\n")

