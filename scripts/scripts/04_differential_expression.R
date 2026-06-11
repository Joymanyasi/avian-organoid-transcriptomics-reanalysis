# ============================================================
# SCRIPT 04: Differential Expression & Kienyeji Gap Analysis
# Project: Avian Organoid Transcriptomics Reanalysis
# Dataset: GSE283090 (Sun et al., 2025, Roslin Institute)
# Author: Joy Manyasi Kabaka | ORCID: 0000-0002-2253-3970
#
# PURPOSE:
# Identify DEGs between broiler and layer organoids, then
# highlight kienyeji-relevant genes to argue that indigenous
# African breeds — adapted to thermal stress and endemic
# pathogens — require their own breed-specific organoid models.
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)

# --- LOAD ANNOTATED OBJECT ---
seurat_obj <- readRDS("data/seurat_annotated.rds")

# Confirm required columns
stopifnot("breed"     %in% colnames(seurat_obj@meta.data))
stopifnot("cell_type" %in% colnames(seurat_obj@meta.data))

# --- PSEUDO-BULK: BROILER VS LAYER ACROSS ALL CELLS ---
Idents(seurat_obj) <- "breed"

deg_all <- FindMarkers(
  seurat_obj,
  ident.1          = "Broiler",
  ident.2          = "Layer",
  min.pct          = 0.25,
  logfc.threshold  = 0.25,
  test.use         = "wilcox"
)

deg_all$gene        <- rownames(deg_all)
deg_all$significant <- deg_all$p_val_adj < 0.05 & abs(deg_all$avg_log2FC) > 0.5

write.csv(deg_all, "results/04_DEG_broiler_vs_layer_all_cells.csv")
cat("Total DEGs (all cells):", sum(deg_all$significant, na.rm = TRUE), "\n")

# --- DEG IN EPITHELIAL CELLS SPECIFICALLY ---
Idents(seurat_obj) <- "cell_type"

epithelial_obj <- subset(
  seurat_obj,
  idents = c("Enterocytes","Goblet_cells","Paneth_cells",
             "Tuft_cells","Enteroendocrine")
)

Idents(epithelial_obj) <- "breed"

deg_epithelial <- FindMarkers(
  epithelial_obj,
  ident.1         = "Broiler",
  ident.2         = "Layer",
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

deg_epithelial$gene        <- rownames(deg_epithelial)
deg_epithelial$significant <- deg_epithelial$p_val_adj < 0.05 &
                              abs(deg_epithelial$avg_log2FC) > 0.5

write.csv(deg_epithelial, "results/04_DEG_broiler_vs_layer_epithelial.csv")
cat("Epithelial DEGs:", sum(deg_epithelial$significant, na.rm = TRUE), "\n")

# --- VOLCANO PLOT ---
volcano <- ggplot(
  deg_all,
  aes(x = avg_log2FC,
      y = -log10(p_val_adj),
      color = significant)
) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "#D62728"),
    labels = c("Not significant", "Significant")
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "black") +
  theme_minimal(base_size = 13) +
  labs(
    title    = "Broiler vs Layer — Differential Gene Expression",
    subtitle = "GSE283090 | Sun et al. 2025 | All cell types",
    x        = "Log2 Fold Change (Broiler / Layer)",
    y        = "-Log10 Adjusted P-value",
    color    = NULL,
    caption  = "Significant = padj < 0.05 and |LFC| > 0.5"
  )

ggsave("results/04_volcano_broiler_vs_layer.png", volcano, width = 10, height = 7)

# --- KIENYEJI-RELEVANT GENE PANEL ---
# Genes important for thermal stress, pathogen response,
# and innate immunity — traits that define kienyeji adaptation
kienyeji_genes <- c(
  "ANP32A",    # Influenza host restriction
  "MX1",       # Interferon-induced antiviral
  "HSP90AA1",  # Heat shock — thermal adaptation
  "HSPA8",     # Heat shock protein
  "TLR4",      # Pathogen recognition receptor
  "TLR7",      # RNA virus sensing
  "IL6",       # Pro-inflammatory cytokine
  "IFNG",      # Interferon gamma
  "IRF7"       # Interferon regulatory factor
)

# Keep only genes present in the dataset
kienyeji_present <- kienyeji_genes[kienyeji_genes %in% rownames(seurat_obj)]
cat("Kienyeji-relevant genes found:", paste(kienyeji_present, collapse = ", "), "\n")

# Subset DEG table to kienyeji genes
kienyeji_deg <- deg_all %>%
  filter(gene %in% kienyeji_present) %>%
  arrange(p_val_adj)

write.csv(kienyeji_deg, "results/04_kienyeji_relevant_DEGs.csv")

# Dot plot of kienyeji-relevant genes by breed
Idents(seurat_obj) <- "breed"

p_dot <- DotPlot(
  seurat_obj,
  features = kienyeji_present
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Kienyeji-Relevant Genes: Broiler vs Layer") +
  labs(caption = "Dot size = % cells expressing | Color = mean expression")

ggsave("results/04_kienyeji_genes_dotplot.png", p_dot, width = 12, height = 5)

# --- TOP 20 DEGs SUMMARY ---
top20 <- deg_all %>%
  filter(significant == TRUE) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)

write.csv(top20, "results/04_top20_DEGs.csv")

cat("Script 04 complete.\n")
cat("All results saved to results/ folder.\n")
cat("Next: interpret findings in docs/findings_summary.md\n")
