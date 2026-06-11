# Results Folder

## Important Note
This folder will be populated after the analysis pipeline 
is executed. Results are not yet uploaded as the dataset 
download and pipeline execution are pending.

---

## Expected Outputs

| File | Script | Description |
|------|--------|-------------|
| 01_QC_before_filtering.png | 01_quality_control.R | QC metrics all 8 samples before filtering |
| 01_QC_after_filtering.png | 01_quality_control.R | QC metrics all 8 samples after filtering |
| 02_elbow_plot.png | 02_clustering_annotation.R | PCA elbow plot |
| 02_UMAP_clusters.png | 02_clustering_annotation.R | Cell clusters UMAP |
| 02_UMAP_breed.png | 02_clustering_annotation.R | Broiler vs Layer UMAP |
| 02_UMAP_timepoint.png | 02_clustering_annotation.R | D0 vs D3 UMAP |
| 02_UMAP_annotated.png | 02_clustering_annotation.R | Annotated cell types UMAP |
| 02_marker_dotplot.png | 02_clustering_annotation.R | Marker gene dotplot |
| 03_ANP32A_all_cells_UMAP.png | 03_ANP32A_analysis.R | ANP32A across all cells |
| 03_ANP32A_broiler_vs_layer.png | 03_ANP32A_analysis.R | ANP32A breed comparison |
| 03_ANP32A_epithelial_only.png | 03_ANP32A_analysis.R | ANP32A in epithelial cells |
| 03_ANP32A_celltype_breed_dotplot.png | 03_ANP32A_analysis.R | ANP32A per cell type by breed |
| 03_ISG_dotplot.png | 03_ANP32A_analysis.R | Interferon-stimulated genes |
| 03_ANP32A_mean_expression_breed_celltype.csv | 03_ANP32A_analysis.R | Mean ANP32A summary table |
| 04_volcano_broiler_vs_layer.png | 04_differential_expression.R | DEG volcano plot |
| 04_DEG_broiler_vs_layer_all_cells.csv | 04_differential_expression.R | Full DEG table all cells |
| 04_DEG_broiler_vs_layer_epithelial.csv | 04_differential_expression.R | DEG table epithelial cells |
| 04_kienyeji_relevant_DEGs.csv | 04_differential_expression.R | Kienyeji-relevant gene DEGs |
| 04_kienyeji_genes_dotplot.png | 04_differential_expression.R | Kienyeji gene panel dotplot |
| 04_top20_DEGs.csv | 04_differential_expression.R | Top 20 DEGs summary |

---

## Status
Results pending — awaiting data download and pipeline execution.
