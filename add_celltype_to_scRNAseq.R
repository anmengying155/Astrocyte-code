# ============================================================================
# Add Cell Type Annotation to scRNA-seq Data
# Purpose: Load 10X data, annotate cell types, and export macrophage/neutrophil subsets
# Input: Cell Ranger output, cell type annotation file
# Output: *_macrophage_neutrophil.rds files
# ============================================================================

library(Seurat)
library(dplyr)
library(readxl)

# Define paths relative to project root
SCRIPT_DIR <- dirname(parent.frame(2)$ofile)
PROJECT_ROOT <- ifelse(is.null(SCRIPT_DIR), ".", dirname(SCRIPT_DIR))

MATRIX_DIR <- file.path(PROJECT_ROOT, "..", "..", "Seq2_3_Comparison", "summary (2)", "summary", "2_Expression_result", "filtered_feature_bc_matrix")
ORIG_IDENT_CSV <- file.path(PROJECT_ROOT, "..", "..", "Seq2_3_Comparison", "summary (2)", "summary", "2_Expression_result", "orig.ident.csv")
CLUSTER_CSV <- file.path(PROJECT_ROOT, "..", "..", "Seq2_3_Comparison", "summary (2)", "summary", "2_Expression_result", "Loupe_Browser_File", "seurat_loupe_cell_umap_cluster.csv")
CELLTYPE_EXCEL <- file.path(PROJECT_ROOT, "..", "Data", "CellType_Annotation.xlsx")

data <- Read10X(MATRIX_DIR)
seurat_obj <- CreateSeuratObject(counts = data, project = "scRNAseq")

orig_df <- read.csv(ORIG_IDENT_CSV) %>% rename(Sample = orig.ident)
cluster_df <- read.csv(CLUSTER_CSV) %>% rename(Barcode = Cell)
celltype_df <- read_excel(CELLTYPE_EXCEL) %>% rename(Cluster = cluster, CellType = cells)

cluster_df$Cluster <- as.character(cluster_df$Cluster)
celltype_df$Cluster <- as.character(celltype_df$Cluster)

meta_data <- data.frame(Barcode = rownames(seurat_obj@meta.data), seurat_obj@meta.data, stringsAsFactors = FALSE)
meta_data <- meta_data %>%
  left_join(orig_df, by = "Barcode") %>%
  left_join(cluster_df, by = "Barcode") %>%
  left_join(celltype_df, by = "Cluster")
rownames(meta_data) <- meta_data$Barcode
meta_data$Barcode <- NULL
seurat_obj@meta.data <- meta_data
seurat_obj$Cluster <- as.character(seurat_obj$Cluster)
seurat_obj$CellType[is.na(seurat_obj$CellType)] <- "Unknown"

OUTPUT_DIR <- file.path(PROJECT_ROOT, "..", "..", "Seq2_3_Comparison", "summary (2)", "summary", "2_Expression_result")
samples <- c("A0_EV", "A1_EV", "A2_EV")
for (sample_name in samples) {
  cat("\n========== Processing:", sample_name, "==========\n")
  seurat_sample <- seurat_obj[, seurat_obj$Sample == sample_name]
  seurat_subset <- seurat_sample[, seurat_sample$CellType %in% c("Macrophage", "Neutrophil")]
  output_file <- file.path(OUTPUT_DIR, paste0(sample_name, "_macrophage_neutrophil.rds"))
  saveRDS(seurat_subset, output_file)
  cat("Cells in", sample_name, "Macrophage+Neutrophil:", ncol(seurat_subset), "\n")
  cat("CellType distribution:\n")
  print(table(seurat_subset$CellType))
}