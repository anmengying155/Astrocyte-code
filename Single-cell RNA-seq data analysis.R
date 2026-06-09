# -------------------------------
# 1. 加载依赖包
# -------------------------------
if (!require("Seurat")) install.packages("Seurat")
if (!require("patchwork")) install.packages("patchwork")
if (!require("cluster")) install.packages("cluster")  # 层次聚类
if (!require("Rtsne")) install.packages("Rtsne")      # t-SNE
library(Seurat)
library(patchwork)
library(cluster)
library(Rtsne)
library(ggplot2)

# -------------------------------
# 2. 创建输出文件夹
# -------------------------------
output_dir <- "A0-A2-scRNA_output"
image_dir <- file.path(output_dir, "images")
data_dir <- file.path(output_dir, "data")

# 检查文件夹是否存在，若不存在则创建
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

if (!dir.exists(image_dir)) {
  dir.create(image_dir)
}

if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
# -------------------------------
# 3. 读取10X数据
# -------------------------------
data_path <- "F:/测序结果/第三次-星胶EV测序/Summary (1)/Summary/2_Expression_result/filtered_feature_bc_matrix" 
source_data <- Read10X(data.dir = data_path)
seurat_raw <- CreateSeuratObject(counts = source_data, project = "10X_Data", min.cells = 3, min.features = 200)
# 对 Seurat 对象中的数据进行标准化
seurat_raw <- NormalizeData(seurat_raw, normalization.method = "LogNormalize", scale.factor = 10000)
# -------------------------------
# 4. 添加线粒体基因比例 & 质控
# -------------------------------
seurat_raw[["percent.mt"]] <- PercentageFeatureSet(seurat_raw, pattern = "^mt-")  # 小鼠用 ^mt-

VlnPlot(seurat_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(image_dir, "quality_control_vlnplot.png"))

seurat_qc <- subset(seurat_raw, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

# -------------------------------
# 5. 标准化 & 筛选高变基因
# -------------------------------
seurat_norm <- NormalizeData(seurat_qc, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_norm <- FindVariableFeatures(seurat_norm, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(seurat_norm), 10)
VariableFeaturePlot(seurat_norm) + LabelPoints(points = top10, repel = TRUE)
ggsave(file.path(image_dir, "variable_feature_plot.png"))

# -------------------------------
# 6. 缩放 & PCA降维
# -------------------------------
seurat_pca <- ScaleData(seurat_norm)
seurat_pca <- RunPCA(seurat_pca, features = VariableFeatures(seurat_norm))

VizDimLoadings(seurat_pca, dims = 1:2, reduction = "pca")
ggsave(file.path(image_dir, "pca_dimplot.png"))

DimPlot(seurat_pca, reduction = "pca")
ggsave(file.path(image_dir, "pca_clusters.png"))

ElbowPlot(seurat_pca)
ggsave(file.path(image_dir, "pca_elbowplot.png"))

# -------------------------------
# 7. 聚类 & UMAP可视化
# -------------------------------
seurat_umap <- FindNeighbors(seurat_pca, dims = 1:10)
seurat_umap <- FindClusters(seurat_umap, resolution = 0.8)
seurat_umap <- RunUMAP(seurat_umap, dims = 1:10)

DimPlot(seurat_umap, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave(file.path(image_dir, "umap_clusters.png"))

# -------------------------------
# 8. t-SNE 可视化
# -------------------------------
seurat_tsne <- RunTSNE(seurat_umap, dims = 1:10)
DimPlot(seurat_tsne, reduction = "tsne", label = TRUE, pt.size = 0.5)
ggsave(file.path(image_dir, "tsne_clusters.png"))

# -------------------------------
# 10. Marker基因识别
# -------------------------------
#markers <- FindAllMarkers(seurat_umap, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.26)
#write.csv(markers, file.path(data_dir, "cluster_markers.csv"))

# 可视化常见免疫细胞 marker
#FeaturePlot(seurat_umap, features = c("Cd3"))
#ggsave(file.path(image_dir, "immune_marker_featureplot.png"))

# -------------------------------
# 11. 保存分析对象
# -------------------------------
saveRDS(seurat_raw, file = file.path(data_dir, "seurat_raw.rds"))
saveRDS(seurat_qc, file = file.path(data_dir, "seurat_qc.rds"))
saveRDS(seurat_norm, file = file.path(data_dir, "seurat_norm.rds"))
saveRDS(seurat_pca, file = file.path(data_dir, "seurat_pca.rds"))
saveRDS(seurat_umap, file = file.path(data_dir, "seurat_umap.rds"))
saveRDS(seurat_tsne, file = file.path(data_dir, "seurat_tsne.rds"))

# -------------------------------
# 完成提示
# -------------------------------
cat("✅ 分析已完成，结果已保存至文件夹：", output_dir, "\n")


