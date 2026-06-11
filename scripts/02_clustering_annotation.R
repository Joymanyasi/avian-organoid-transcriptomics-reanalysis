# ============================================================
# SCRIPT 02: Normalization, Clustering & Cell Type Annotation
# Project: Avian Organoid Transcriptomics Reanalysis
# Dataset: GSE283090 (Sun et al., 2025, Roslin Institute)
# Author: Joy Manyasi Kabaka | ORCID: 0000-0002-2253-3970
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# --- LOAD FILTERED OBJECT ---
seurat_obj <- readRDS("data/seurat_merged_filtered.rds")

# --- NORMALIZE ---
seurat_obj <- NormalizeData(seurat_obj)

# --- VARIABLE FEATURES ---
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures        = 2000
)

# --- SCALE ---
seurat_obj <- ScaleData(seurat_obj)

# --- PCA ---
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# Elbow plot to choose optimal number of PCs
p_elbow <- ElbowPlot(seurat_obj, ndims = 30)
ggsave("results/02_elbow_plot.png", p_elbow, width = 7, height = 5)

# --- CLUSTERING ---
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj,  resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj,       dims = 1:20)

# --- UMAP PLOTS ---
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
      ggtitle("Seurat Clusters")

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "breed") +
      ggtitle("Broiler vs Layer")

p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "timepoint") +
      ggtitle("D0 vs D3")

ggsave("results/02_UMAP_clusters.png",   p1, width = 8, height = 6)
ggsave("results/02_UMAP_breed.png",      p2, width = 8, height = 6)
ggsave("results/02_UMAP_timepoint.png",  p3, width = 8, height = 6)

# --- MARKER GENES FROM Sun et al. 2025 ---
# Used to assign cell type identity to each cluster
marker_genes <- c(
  "FABP2",   # Enterocytes
  "MUC2",    # Goblet cells
  "DEFA",    # Paneth cells
  "TRPM5",   # Tuft cells
  "COL1A1",  # Fibroblasts
  "ACTA2",   # Smooth muscle
  "CSF1R",   # Macrophages
  "CD3D"     # T cells
)

p_markers <- FeaturePlot(
  seurat_obj,
  features  = marker_genes,
  ncol      = 4,
  reduction = "umap"
)
ggsave("results/02_marker_genes_UMAP.png", p_markers, width = 16, height = 8)

p_dot <- DotPlot(seurat_obj, features = marker_genes) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
         ggtitle("Marker Gene Expression per Cluster")
ggsave("results/02_marker_dotplot.png", p_dot, width = 12, height = 6)

# --- ASSIGN CELL TYPE LABELS ---
# Cluster numbers will vary with your data.
# After inspecting the UMAP and dot plot above, replace the
# right-hand values below with the cell types you observe.
# Example based on Sun et al. 2025 expected clusters:

cluster_labels <- c(
  "0"  = "Enterocytes",
  "1"  = "Fibroblasts",
  "2"  = "Goblet_cells",
  "3"  = "Smooth_muscle",
  "4"  = "Enterocytes",
  "5"  = "Macrophages",
  "6"  = "Paneth_cells",
  "7"  = "Tuft_cells",
  "8"  = "T_cells",
  "9"  = "Fibroblasts",
  "10" = "Enteroendocrine",
  "11" = "Unknown"
)

seurat_obj <- RenameIdents(seurat_obj, cluster_labels)
seurat_obj$cell_type <- as.character(Idents(seurat_obj))

# UMAP with cell type labels
p_annotated <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("Annotated Cell Types — GSE283090")

ggsave("results/02_UMAP_annotated.png", p_annotated, width = 10, height = 7)

# --- SAVE ---
saveRDS(seurat_obj, "data/seurat_annotated.rds")
cat("Script 02 complete. Cell type column 'cell_type' added and saved.\n")
