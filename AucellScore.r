
# -------------------------------
# 2. 创建输出文件夹
# -------------------------------
output_dir <- "A0-A2_result"
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
# 加载数据
#-------------------------------
# 5. 读取你自己的表达矩阵文件
#-------------------------------
library(Seurat)
library(openxlsx)
# 读取 filtered_feature_bc_matrix 文件夹
load("singleR注释库/sub_Macrophage+Mono_Filter_Basicinfo.rData")
ls()
seurat_obj <-pbmc_filt
exprMatrix <- GetAssayData(seurat_obj, slot = "data")
# 查看表达矩阵信息
dim(exprMatrix)  # 查看维度
head(rownames(exprMatrix))  # 查看基因名
head(colnames(exprMatrix))  # 查看细胞条形码
exprMatrix[1:5, 1:5]  # 查看前几行
seurat_obj<- readRDS("A0-A2_seurat_obj.rds")
# 查看基本信息
# 读入外部计算好的 UMAP 坐标
external_umap <- read.csv("F:/测序结果/第三次-星胶EV测序/Summary (1)/Summary/2_Expression_result/Loupe_Browser_File/seurat_loupe_cell_umap_pos.csv", row.names = 1)
# 确保行名和 seurat_obj 的细胞名一致
head(rownames(external_umap))
head(colnames(seurat_obj))  # 对比一下是否一致
# 把 UMAP 数据变成 matrix
umap_matrix <- as.matrix(external_umap)
head(umap_matrix)
# 加入 Seurat 对象作为自定义降维结果
seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = umap_matrix,
                                                      key = "UMAP_",
                                                      assay = DefaultAssay(seurat_obj))

# -------------------------------
# UMAP可视化
# -------------------------------
DimPlot(seurat_obj, reduction = "umap")
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 3) + 
  ggtitle("seurat_umap") + coord_fixed()
# -------------------------------
# singleR注释
# -------------------------------
MR.se <- MouseRNAseqData()
IGD.se <- ImmGenData()
SingleR_data <- GetAssayData(seurat_obj, slot = "data")
# 对 MR 数据集进行注释
SingleR_data_MR <- SingleR(test = SingleR_data, ref = MR.se, labels = MR.se$label.main)
seurat_obj$MR_labels <- SingleR_data_MR$labels
saveRDS(seurat_obj, file = "A0-A2labled.rds")
# 计算每种细胞类型的细胞数量
cell_counts_MR <- table(seurat_obj$MR_labels)
# -------------------------------
# UMAP 可视化 MR 注释
# -------------------------------
dimplot_MR <- DimPlot(seurat_obj, reduction = "umap", group.by = "MR_labels", label = FALSE, pt.size = 2) + 
  ggtitle("Cell Types - MR Labels") + coord_fixed() +
  theme(text = element_text(size = 12))  # 调整字体大小
#-------------------------------
# 6. 准备基因集（Gene Set）
#-------------------------------
library(AUCell)
library(GSEABase)
# 加载AUCell内置的例子基因集文件（⚠️你也可以换成自己的gmt文件）
#gmtFile <- file.path(system.file('examples', package='AUCell'), "geneSignatures.gmt")
gmtFile<-file.path('D:/File/document/单细胞测序代码/AUcell/referencegeneset.gmt')
# 载入基因集
geneSets <- getGmt(gmtFile)
# 筛选和表达矩阵有交集的基因集
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix))
cbind(nGenes(geneSets))  # 查看每个基因集的基因数
head(geneSets)
#-------------------------------
# 7.计算得分
#-------------------------------
library(parallel)
options(mc.cores = 4)  # 设置使用的核心数

# 获取 orig.ident 信息
exprMatrix <- GetAssayData(seurat_obj, slot = "data")
cells_rankings <- AUCell_buildRankings(exprMatrix)
aucMaxRank<- 1026
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = aucMaxRank)
save(cells_AUC, file = "A0-A2_AUC.RData")

## 案例
library(Seurat)
library(AUCell)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)

# Step 1: Get AUCell matrix and convert to a data frame
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))
auc_df$labels <- seurat_obj$MR_labels[match(rownames(auc_df), colnames(seurat_obj))]
auc_df$orig.ident <- seurat_obj$orig.ident[match(rownames(auc_df), colnames(seurat_obj))]

# Step 3: Reshape the data into a long format
auc_long <- pivot_longer(auc_df, cols = -c(orig.ident,labels),
                         names_to = "gene_set", values_to = "auc_score")

# Step 4: Calculate the mean AUCell for each combination of celltype, and gene_set
mean_df <- auc_long %>%
  group_by(labels, gene_set,orig.ident) %>%
  summarise(mean_auc = mean(auc_score, na.rm = TRUE), .groups = "drop")

# Step 5: Get unique gene sets and sample identifiers
unique_gene_sets <- unique(auc_long$gene_set)
selected_gene_set_indices <- c(6)
selected_gene_sets <- unique_gene_sets[selected_gene_set_indices]
unique_celltypes <- unique(auc_long$labels)
selected_celltypes_indices <- c(3,5,6,8,11,12,13)
selected_celltypes<- unique_celltypes[selected_celltypes_indices]
unique_orig.ident <- unique(auc_long$orig.ident)
selected_orig.ident_indices <- c(2)
selected_orig.ident<- unique_orig.ident[selected_orig.ident_indices]
# Step 6: Filter the data to select the specified gene sets and orig.ident
auc_long_selected <- auc_long %>% filter(gene_set %in% selected_gene_sets & labels %in% selected_celltypes & orig.ident %in% selected_orig.ident )
mean_df_selected <- mean_df%>% filter(gene_set %in% selected_gene_sets & labels %in% selected_celltypes & orig.ident %in% selected_orig.ident )
# Create a data frame for the mean difference
write.csv(as.data.frame(mean_df_selected), "mean_df_selected.csv")

# 根据中位数排序
label_order <- auc_long_selected %>%
  group_by(labels) %>%
  summarize(median_auc = median(auc_score)) %>%
  arrange(desc(median_auc)) %>%
  pull(labels)

# 应用排序
auc_long_selected$labels <- factor(auc_long_selected$labels, levels = label_order)

color_map <- c("MO" = "#FF6347", "REP" = "#4682B4", "NT" = "#32CD32")

    p <- ggplot(auc_long_selected,
              aes(x = auc_score, y = labels , fill = orig.ident, color = orig.ident)) +
    geom_density_ridges(alpha = 0.4, scale = 1.0, rel_min_height = 0.01, size = 0.6) +  # 透明度调整为0.4，外围线条大小为0.6
    theme_minimal(base_size = 13) +
    # 使用定义的颜色映射
    scale_fill_manual(values = color_map) +
    scale_color_manual(values = color_map) +
    labs(title =  "positive_regulation_of_angiogenesis",
         x = "AUCell Score", y = "orig.ident") +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none"
    )
  
  print(p)  # 显示图
  