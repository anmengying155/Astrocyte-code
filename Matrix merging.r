# -------------------------------
# 加载必要的库
# -------------------------------
library(Seurat)
library(SingleR)
library(celldex)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(loupeR)
library(InterCellDB)  # 0.9.2
library(Seurat)       # 3.2.0
library(dplyr)        # 1.0.5
library(future)       # 1.21.0
library(future.apply)
library(openxlsx)
# -------------------------------
# 加载 Seurat 对象
# -------------------------------
setwd("D:/File/document/单细胞测序代码") 
load("D:/File/document/单细胞测序代码/singleR注释库/sub_巨噬+单核_cells_去批次_bimod_Filter_Basicinfo.rData")
seurat_obj <- readRDS("scRNA_output/data/seurat_umap.rds")
unique(pbmc_filt$seurat_clusters)
head(pbmc_filt)
celltype_info <- read.xlsx("D:/File/document/单细胞测序代码/cell_annot (6-Macro)-2.xlsx")
head(celltype_info)
pbmc_filt$celltype <- celltype_info$cells[match(pbmc_filt$seurat_clusters,celltype_info$cluster )]

# -------------------------------
# 选择目标细胞群
# -------------------------------
brain <- subset(pbmc_filt, subset = orig.ident == "A2_EV" & celltype %in% c("MF1", "MF2"))
brain_MF2 <- subset(pbmc_filt, subset = orig.ident == "A2_EV" & celltype %in% c("MF2"))
blood <- subset(seurat_obj, subset = seurat_clusters %in% c("9", "11"))
blood_M<- subset(seurat_obj, subset = seurat_clusters %in% c("5"))
blood$orig.ident <- "Blood"
dim(blood_M)
dim(brain)

# -------------------------------
# 标准化和寻找特征基因
# -------------------------------
brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain)
head(brain)
blood <- NormalizeData(blood)
blood <- FindVariableFeatures(blood)
blood$celltype<- "monocyte"
head(blood)
# -------------------------------
# 整合准备
# -------------------------------
# 创建一个 list
seurat_list <- list(brain, blood)
# 查找锚点
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30)
# 整合数据
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
seurat_integrated <- NormalizeData(seurat_integrated)
# -------------------------------
# 后续分析（降维 + 聚类）
# -------------------------------
# 设置默认数据
DefaultAssay(seurat_integrated) <- "integrated"

# 归一化、主成分分析
seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30)

# 可视化
ElbowPlot(seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:20)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.8)
seurat_integrated$orig.ident <- seurat_integrated$celltype
# 绘图查看整合效果
umap_plot_1=DimPlot(seurat_integrated, reduction = "umap", group.by = "orig.ident", pt.size = 1.5) + ggtitle("Integrated UMAP")
# 绘图查看整合效果
umap_plot_1=DimPlot(seurat_integrated, reduction = "umap", group.by = "celltype", pt.size = 1.5) + ggtitle("Integrated UMAP")
# 使用 ggsave 保存为 SVG 格式
# 使用 ggsave 保存为 SVG 格式
ggsave("Integrated_UMAP4.svg", plot = umap_plot_1, device = "svg", width = 8, height = 6)
umap_plot_2=DimPlot(seurat_integrated, reduction = "umap", label = TRUE, pt.size = 1.5) + ggtitle("Clusters")
ggsave("Clusters_UMAP4.svg", plot = umap_plot_2, device = "svg", width = 8, height = 6)

## ------------------------------------------------------------------------------

##-------------------------------------------------------------------------------
# 可视化这些基因的表达 —— 用小提琴图（分组显示）
VlnPlot(
  seurat_integrated,
  features = top_gene_names,
  group.by = "celltype", # 按照orig.ident分组，如Brain和Blood
  pt.size = 0,             # 不显示散点，更清晰
  ncol = 5                 # 每行5张图
)

DoHeatmap(
  seurat_integrated,
  features = top_gene_names,
  group.by = "celltype" # 分组显示，比如 Brain 和 Blood
) + 
  scale_fill_gradientn(colors = c("#242457", "white", "#8e0505"))  # 可选：设置颜色渐变，更直观


# 提取 normalized 表达矩阵
heatmap_data <- GetAssayData(seurat_integrated, slot = "scale.data")[top20_genes, ]
# 查看前几行
head(heatmap_data)
# 可选：按分组排序
cell_order <- seurat_integrated@meta.data %>%
  dplyr::arrange(orig.ident) %>%
  rownames()
heatmap_data <- heatmap_data[, cell_order]
# 保存为 CSV 文件
write.csv(as.data.frame(heatmap_data), "Heatmap_Expression_Top20 MF1 and MF2_vs_Blood_scaled.csv")



library(monocle3)
library(SeuratWrappers)
# 假设你已经合并好 Seurat 对象为 combined
# 并已完成标准的 NormalizeData、FindVariableFeatures、RunPCA、RunUMAP 等流程

# 转换为 CellDataSet
cds <- as.cell_data_set(seurat_integrated)

# 将原始分群和来源样本信息加进去
cds@colData$cluster <- seurat_integrated$celltype
cds@colData$sample <- seurat_integrated$orig.ident  # 或其他能区分 blood / brain 的字段

# 如果之前已经降维过（如 PCA），这一步可以直接继续
cds <- preprocess_cds(cds, num_dim = 30)

# 降维和聚类（用UMAP）
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)

# 构建轨迹图
cds <- learn_graph(cds)
# 查看 cluster 和 sample 的分布，选择 blood 中的 cluster 作为起点
table(colData(cds)$cluster, colData(cds)$sample)

# 假设 blood 中的 cluster "1" 是分化起点
start_cell <- colnames(cds)[which(colData(cds)$cluster == "1" & colData(cds)$sample == "blood")[1]]

# 设置 pseudotime 起始细胞
cds <- order_cells(cds, root_cells = start_cell)
# 按 pseudotime 上色
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)

# 按 sample 来源上色（例如 blood vs brain）
plot_cells(
  cds,
  color_cells_by = "sample"
)

# 按 cluster 上色
plot_cells(
  cds,
  color_cells_by = "cluster"
)
# 示例：选择3个基因展示它们在 pseudotime 上的表达趋势
plot_genes_in_pseudotime(cds[c("GeneA", "GeneB", "GeneC"), ])
