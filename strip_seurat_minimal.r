strip_seurat_minimal <- function(seurat_obj) {
  # 只保留 RNA assay 和 counts 数据
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj@assays <- seurat_obj@assays["RNA"]
  seurat_obj@reductions <- list()
  seurat_obj@graphs <- list()
  seurat_obj@neighbors <- list()
  seurat_obj@commands <- list()
  
  # 只保留 orig.ident 信息
  if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data <- seurat_obj@meta.data[, "orig.ident", drop = FALSE]
  } else {
    seurat_obj@meta.data <- data.frame(orig.ident = rep("unknown", ncol(seurat_obj)))
  }
  
  return(seurat_obj)
}