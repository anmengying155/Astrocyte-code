# -------------------------------
# 加载必要的库
# -------------------------------
library(SingleR)
library(celldex)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(dplyr)

# -------------------------------
# 2. 创建输出文件夹
# -------------------------------
output_dir <- "sigleR_result"
image_dir <- file.path(output_dir, "images")
data_dir <- file.path(output_dir, "data")

if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(image_dir)) dir.create(image_dir)
if (!dir.exists(data_dir)) dir.create(data_dir)

# -------------------------------
# 加载数据
# -------------------------------
seurat_obj <- readRDS("scRNA_output/data/seurat_umap.rds")

# -------------------------------
# 获取参考数据集
# -------------------------------
MR.se <- MouseRNAseqData()
IGD.se <- ImmGenData()

# -------------------------------
# singleR 注释
# -------------------------------
SingleR_data <- GetAssayData(seurat_obj, slot = "data")

SingleR_data_MR <- SingleR(test = SingleR_data, ref = MR.se, labels = MR.se$label.fine)
SingleR_data_IGD <- SingleR(test = SingleR_data, ref = IGD.se, labels = IGD.se$label.fine)

# -------------------------------
# 添加注释标签到 Seurat 对象
# -------------------------------
seurat_obj$MR_labels <- SingleR_data_MR$labels
seurat_obj$IGD_labels <- SingleR_data_IGD$labels

# -------------------------------
# UMAP 可视化 MR 和 IGD 注释
# -------------------------------
dimplot_MR <- DimPlot(seurat_obj, reduction = "umap", group.by = "MR_labels", label = FALSE, pt.size = 2) + 
  ggtitle("Cell Types - MR Labels") + coord_fixed() + theme(text = element_text(size = 12))
ggsave(file.path(image_dir, "dimplot_MR.png"), plot = dimplot_MR, width = 10, height = 8, dpi = 300)

dimplot_IGD <- DimPlot(seurat_obj, reduction = "umap", group.by = "IGD_labels", label = FALSE, pt.size = 2) + 
  ggtitle("Cell Types - IGD Labels") + coord_fixed() + theme(text = element_text(size = 12))
ggsave(file.path(image_dir, "dimplot_IGD.png"), plot = dimplot_IGD, width = 10, height = 8, dpi = 300)

# -------------------------------
# Marker基因提取函数
# -------------------------------
extract_markers <- function(seurat_obj, label_col, output_name) {
  Idents(seurat_obj) <- seurat_obj[[label_col]]
  markers <- FindAllMarkers(seurat_obj, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers, file.path(data_dir, paste0(output_name, "_marker_genes.csv")), row.names = FALSE)
  return(markers)
}

# -------------------------------
# 提取 MR_labels Marker
# -------------------------------
markers_MR <- extract_markers(seurat_obj, "MR_labels", "MR_labels")

# -------------------------------
# 提取 IGD_labels Marker
# -------------------------------
markers_IGD <- extract_markers(seurat_obj, "IGD_labels", "IGD_labels")

# -------------------------------
# 示例：绘制MR_labels的top10 marker基因热图
# -------------------------------
top10_MR <- markers_MR %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>% 
  pull(gene) %>% 
  unique()

# 提取表达矩阵
expr_mat <- GetAssayData(seurat_obj, slot = "data")[top10_MR, ]

# 对细胞分组排序
cell_order <- order(seurat_obj$MR_labels)
annotation_col <- data.frame(CellType = seurat_obj$MR_labels)
rownames(annotation_col) <- colnames(seurat_obj)

# 绘图
pheatmap(expr_mat[, cell_order],
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         scale = "row",
         fontsize = 8,
         main = "Top10 Marker Genes by MR_labels")

# -------------------------------
# 结束
# -------------------------------
