library(Seurat)
RDS_FILE <- r"(d:\单细胞测序\Seq2_3_Comparison\summary (2)\summary\2_Expression_result\filtered_feature_bc_matrix\macrophage_neutrophil.rds)"
seurat_obj <- readRDS(RDS_FILE)
cat("=== orig.ident values ===\n")
print(unique(seurat_obj@meta.data$orig.ident))
cat("\n=== Sample distribution ===\n")
print(table(seurat_obj@meta.data$orig.ident))
cat("\n=== Barcode examples ===\n")
print(head(rownames(seurat_obj@meta.data)))