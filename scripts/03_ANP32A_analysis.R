# ============================================================
# SCRIPT 03: ANP32A Host-Restriction Locus Analysis
# Project: Avian Organoid Transcriptomics Reanalysis
# Dataset: GSE283090 (Sun et al., 2025, Roslin Institute)
# Author: Joy Manyasi Kabaka | ORCID: 0000-0002-2253-3970
#
# SCIENTIFIC RATIONALE:
# ANP32A is a critical host-restriction factor for avian influenza.
# Variant p.N129D determines susceptibility to LPAIV/NDV.
# The White Leghorn reference genome carries one allele — masking
# allelic diversity in non-commercial breeds like kienyeji.
# This script asks: does ANP32A expression already differ between
# commercial broiler and layer breeds? If yes, the gap to kienyeji
# is even greater and breed-specific models are essential.
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# --- LOAD ANNOTATED OBJECT ---
# Requires cell_type column created in script 02
seurat_obj <- readRDS("data/seurat_annotated.rds")

# Confirm cell_type column exists
stopifnot("cell_type" %in% colnames(seurat_obj@meta.data))

# --- ANP32A EXPRESSION ACROSS ALL CELLS ---
p1 <- FeaturePlot(
  seurat_obj,
  features  = "ANP32A",
  reduction = "umap"
) + ggtitle("ANP32A Expression — All Cell Types")

ggsave("results/03_ANP32A_all_cells_UMAP.png", p1, width = 8, height = 6)

# --- ANP32A BY BREED ---
p2 <- VlnPlot(
  seurat_obj,
  features = "ANP32A",
  group.by = "breed",
  pt.size  = 0
) + ggtitle("ANP32A Expression: Broiler vs Layer — All Cells")

ggsave("results/03_ANP32A_broiler_vs_layer.png", p2, width = 7, height = 5)

# --- ANP32A IN EPITHELIAL CELLS ONLY ---
# Epithelial cells are primary entry point for respiratory pathogens
epithelial_cells <- subset(
  seurat_obj,
  subset = cell_type %in% c("Enterocytes", "Goblet_cells",
                             "Paneth_cells", "Tuft_cells",
                             "Enteroendocrine")
)

p3 <- VlnPlot(
  epithelial_cells,
  features = "ANP32A",
  group.by = "breed",
  pt.size  = 0
) + ggtitle("ANP32A in Epithelial Cells: Broiler vs Layer")

ggsave("results/03_ANP32A_epithelial_only.png", p3, width = 7, height = 5)

# --- ANP32A ACROSS ALL CELL TYPES BY BREED ---
p4 <- DotPlot(
  seurat_obj,
  features = "ANP32A",
  group.by = "cell_type",
  split.by = "breed",
  cols     = c("blue", "red")
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("ANP32A per Cell Type: Broiler (blue) vs Layer (red)")

ggsave("results/03_ANP32A_celltype_breed_dotplot.png", p4, width = 12, height = 5)

# --- INTERFERON STIMULATED GENES ---
isg_genes <- c("MX1","IFIT1","IFIT3","ISG15","OAS1","IRF7")

# Check which ISGs are present in this dataset
isg_present <- isg_genes[isg_genes %in% rownames(seurat_obj)]
cat("ISGs found in dataset:", paste(isg_present, collapse = ", "), "\n")

if (length(isg_present) > 0) {
  p5 <- DotPlot(
    seurat_obj,
    features = isg_present,
    group.by = "breed"
  ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Interferon-Stimulated Genes: Broiler vs Layer")
  
  ggsave("results/03_ISG_dotplot.png", p5, width = 10, height = 5)
}

# --- SUMMARY TABLE: MEAN ANP32A PER BREED PER CELL TYPE ---
avg_expr <- AverageExpression(
  seurat_obj,
  features = "ANP32A",
  group.by = c("cell_type", "breed"),
  assays   = "RNA"
)

write.csv(
  as.data.frame(avg_expr$RNA),
  "results/03_ANP32A_mean_expression_breed_celltype.csv"
)

cat("Script 03 complete.\n")
cat("Key question answered: Does ANP32A differ between broiler and layer?\n")
cat("Check results/03_ANP32A_broiler_vs_layer.png\n")
