
# -------------------------------
# 2. 创建输出文件夹
# -------------------------------
output_dir <- "Tcell_result"
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
load("D:/File/document/单细胞测序代码/singleR注释库/sub_T_cells_去批次_bimod_Filter_Basicinfo.rData")
ls()
seurat_obj <-pbmc_filt
exprMatrix <- GetAssayData(seurat_obj, slot = "data")
# 查看表达矩阵信息
dim(exprMatrix)  # 查看维度
head(rownames(exprMatrix))  # 查看基因名
head(colnames(exprMatrix))  # 查看细胞条形码
exprMatrix[1:5, 1:5]  # 查看前几行

# 查看基本信息
head(seurat_obj)
unique(seurat_obj$seurat_clusters)
unique(seurat_obj$celltype)
cluster_info <- read.xlsx("F:/测序结果/第2和3对比/sub_T_cells_去批次_bimod/cell_clusterinfo (6-Tcell).xlsx")
head(cluster_info)
# 将 cluster_info 中的细胞类型赋值到 seurat_obj$celltype 中
seurat_obj$celltype <- cluster_info$cells[match(seurat_obj$seurat_clusters, cluster_info$cluster)]
head(seurat_obj)
# 检查是否有 NA 值
table(is.na(seurat_obj$celltype))
# 取出非 NA 的细胞名
valid_cells <- colnames(seurat_obj)[!is.na(seurat_obj$celltype)]
# 用细胞名子集化 Seurat 对象
seurat_obj <- subset(seurat_obj, cells = valid_cells)
# 查看 Seurat 对象的所有元数据列
colnames(seurat_obj@meta.data)
# 4. 查看添加后的 Seurat 对象的元数据
head(seurat_obj@meta.data)
dim(seurat_obj)
unique(seurat_obj$orig.ident)
head(seurat_obj)
seurat_obj <- NormalizeData(seurat_obj) 
exprMatrix <- GetAssayData(seurat_obj, slot = "data")
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
A2_EV_obj <- subset(seurat_obj, orig.ident == "A2_EV")
exprMatrix <- GetAssayData(A2_EV_obj, slot = "data")
cells_rankings <- AUCell_buildRankings(exprMatrix)
aucMaxRank<- 803
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = aucMaxRank)
save(cells_AUC, file = paste0(data_dir, "A2_EV_AUC.RData"))

# 获取 orig.ident 信息
DBCO_obj <- subset(seurat_obj, orig.ident == "DBCO")
exprMatrix <- GetAssayData(DBCO_obj, slot = "data")
cells_rankings <- AUCell_buildRankings(exprMatrix)
aucMaxRank<- 1754
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = aucMaxRank)
save(cells_AUC, file = paste0(data_dir, "DBCO_AUC.RData"))

## 案例
library(Seurat)
library(AUCell)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)

load("Tcell_result/dataA2_EV_AUC.RData")
# Step 1: Get AUCell matrix and convert to a data frame
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))
auc_df$celltype <- A2_EV_obj $celltype[match(rownames(auc_df), colnames(A2_EV_obj))]


load("Tcell_result/dataDBCO_AUC.RData")
# Step 1: Get AUCell matrix and convert to a data frame
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))
auc_df$celltype <- DBCO_obj $celltype[match(rownames(auc_df), colnames(DBCO_obj))]



# Step 3: Reshape the data into a long format
auc_long <- pivot_longer(auc_df, cols = -c(celltype),
                         names_to = "gene_set", values_to = "auc_score")

# Step 4: Calculate the mean AUCell for each combination of celltype, and gene_set
mean_df <- auc_long %>%
  group_by(celltype, gene_set) %>%
  summarise(mean_auc = mean(auc_score, na.rm = TRUE), .groups = "drop")

# Step 5: Get unique gene sets and sample identifiers
unique_gene_sets <- unique(mean_df$gene_set)
selected_gene_set_indices <- c(2, 3, 4, 6, 7, 9, 12)
selected_gene_sets <- unique_gene_sets[selected_gene_set_indices]
unique_celltypes <- unique(mean_df$celltype)
# Step 6: Filter the data to select the specified gene sets and orig.ident
mean_df_selected <- mean_df %>% filter(gene_set %in% selected_gene_sets )
auc_long_selected <- auc_long %>% filter(gene_set %in% selected_gene_sets)
# Step 9: Calculate the mean differences between experimental groups and DBCO (control group)
# Create a data frame for the mean difference

# 2. 数据转换为宽格式
heatmap_data <- mean_df_selected %>%
  pivot_wider(names_from = celltype, values_from = mean_auc)

# 3. 将 gene_set 设置为行名，并转换为 data.frame
heatmap_data <- as.data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data$gene_set
heatmap_data <- heatmap_data %>% dplyr::select(-gene_set)

# 4. 标准化数据（行或列）#00ff26#cbf7d1#cbf7e8##1a1a53ea#c4ffe6#bceed9
heatmap_scaled_data <- t(scale(t(heatmap_data)))  # 按行标准化（每个基因）

png("heatmap_output.png", width = 1000, height = 800)  # 设置文件名和图形的宽高

# 绘制热图，使用从 royalblue 到 firebrick3 的颜色渐变
pheatmap(heatmap_data_scaled, 
         scale = "none",  # 如果已经标准化，则不需要再次标准化
         clustering_distance_rows = "euclidean",  # 行的聚类距离
         clustering_distance_cols = "euclidean",  # 列的聚类距离
         color = colorRampPalette(c("royalblue", "white", "firebrick3"))(56),
         show_rownames = TRUE,  # 显示行名（gene_set）
         show_colnames = TRUE,  # 显示列名（celltype）
         cluster_rows = FALSE,  # 不对行进行聚类
         cluster_cols = FALSE,  # 不对列进行聚类
         cellwidth = 100,  # 设置单元格的宽度
         cellheight = 100)  # 设置单元格的高度，使方格正方形

# 关闭设备并保存图像
dev.off()
write.csv(heatmap_data, "DBCO_heatmap_data.csv", row.names = TRUE)
write.csv(heatmap_scaled_data, "DBCO_heatmap_scaled_data.csv", row.names = TRUE)






