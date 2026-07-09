#-------------------------------
# 1. 安装需要的 R/Bioconductor 包
#-------------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("AUCell")
# 并行计算支持
BiocManager::install(c("doMC", "doRNG", "doSNOW"))
# 主分析所需包
BiocManager::install(c("mixtools", "GEOquery", "SummarizedExperiment"))
# 后续可视化等扩展包（可选）
BiocManager::install(c("DT", "plotly", "NMF", "d3heatmap", "shiny", "rbokeh",
                       "dynamicTreeCut", "R2HTML", "Rtsne", "zoo"))
browseVignettes("AUCell")
#-------------------------------
# 2. 设置工作目录（可以自定义路径）
#-------------------------------
dir.create("AUCell_2")        # 你可以更换目录名
setwd("AUCell_2")
#-------------------------------
# 3. 下载 GEO 数据：GSE60361（小鼠单细胞脑表达数据）
#-------------------------------
install.packages("tzdb")
library(GEOquery)
# 可能因为网络不稳定会失败，这里尝试最多20次
attemptsLeft <- 20
while(attemptsLeft > 0) {
  geoFile <- tryCatch(getGEOSuppFiles("GSE60361", makeDirectory=FALSE), error=identity)
  if(methods::is(geoFile, "error")) {
    attemptsLeft <- attemptsLeft - 1
    Sys.sleep(5)
  } else {
    attemptsLeft <- 0
  }
}
#-------------------------------
# 4. 解压和读取表达数据
#-------------------------------
# 自动识别压缩文件名
gzFile <- grep(".txt.gz", basename(rownames(geoFile)), fixed=TRUE, value=TRUE)
txtFile <- gsub(".gz", "", gzFile, fixed=TRUE)
# 解压数据文件
R.utils::gunzip(filename=gzFile, destname=txtFile, remove=TRUE)
# 读取表达矩阵
library(data.table)
geoData <- fread(txtFile, sep="\t")
geneNames <- unname(unlist(geoData[,1, with=FALSE]))
exprMatrix <- as.matrix(geoData[,-1, with=FALSE])
rownames(exprMatrix) <- geneNames
rm(geoData)
# 查看维度
dim(exprMatrix)
exprMatrix[1:5, 1:4]
# 删除原始txt文件节省空间
file.remove(txtFile)
# 保存完整表达矩阵用于之后复用
mouseBrainExprMatrix <- exprMatrix
save(mouseBrainExprMatrix, file="exprMatrix_AUCellVignette_MouseBrain.RData")
#抽样5000个基因用于演示（可自定义数量）
set.seed(333)
exprMatrix <- mouseBrainExprMatrix[sample(rownames(mouseBrainExprMatrix), 5000),]
dim(exprMatrix)


#-------------------------------
# 5. 读取你自己的表达矩阵文件
#-------------------------------
install.packages("parallelly")
library(Seurat)
# 读取 filtered_feature_bc_matrix 文件夹
load("D:/File/document/单细胞测序代码/singleR注释库/sub_T_cells_去批次_bimod_Filter_Basicinfo.rData")
data_dir <- "F:/测序结果/第2和3对比/summary (2)/summary/2_Expression_result/filtered_feature_bc_matrix"
exprMatrix  <- Read10X(data.dir = data_dir)
# 查看表达矩阵信息
dim(exprMatrix)  # 查看维度
head(rownames(exprMatrix))  # 查看基因名
head(colnames(exprMatrix))  # 查看细胞条形码
exprMatrix[1:5, 1:5]  # 查看前几行
seurat_obj <- CreateSeuratObject(counts = exprMatrix)
# 查看基本信息
head(seurat_obj)
barcode_info <- read.csv("F:/测序结果/第2和3对比/summary (2)/summary/2_Expression_result/barcode_orig.ident.csv")
barcode_info$Barcode <- gsub(" ", "", barcode_info$Barcode)
seurat_obj$orig.ident <- barcode_info$orig.ident[match(Cells(seurat_obj), barcode_info$Barcode)]
barcode_info <- read.csv("F:/测序结果/第2和3对比/summary (2)/summary/2_Expression_result/barcode_celltype.csv")
barcode_info$Barcode <- gsub(" ", "", barcode_info$Barcode)
seurat_obj$celltype <- barcode_info$celltype[match(Cells(seurat_obj), barcode_info$Barcode)]
# 取出非 NA 的细胞名
valid_cells <- colnames(seurat_obj)[!is.na(seurat_obj$celltype)]
# 用细胞名子集化 Seurat 对象
seurat_obj <- subset(seurat_obj, cells = valid_cells)
# 查看 Seurat 对象的所有元数据列
colnames(seurat_obj@meta.data)
# 4. 查看添加后的 Seurat 对象的元数据
head(seurat_obj@meta.data)
table(is.na(seurat_obj$orig.ident))  # 看看多少 NA
table(is.na(seurat_obj$celltype))  # 看看多少 NA
dim(barcode_info)
dim(seurat_obj)
unique(seurat_obj$orig.ident)
seurat_obj <- NormalizeData(seurat_obj) 
# 保存 Seurat 对象
saveRDS(seurat_obj, file = "all_seurat_obj.rds")
seurat_obj <- readRDS("D:/File/document/单细胞测序代码/AUCell_tutorial/all_seurat_obj.rds") 
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
# 获取表达数据
seurat_obj <- readRDS("D:/File/document/单细胞测序代码/AUCell_tutorial/all_seurat_obj.rds") 
exprMatrix <- GetAssayData(seurat_obj, slot = "data")
# 获取 orig.ident 信息
A2_EV_obj <- subset(seurat_obj, orig.ident == "A2_EV")
exprMatrix <- GetAssayData(A2_EV_obj, slot = "data")
cells_rankings <- AUCell_buildRankings(exprMatrix)
aucMaxRank<- 1000
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = aucMaxRank)
save(cells_AUC, file="A2_EV_AUC.RData")

# 获取 orig.ident 信息
DBCO_obj <- subset(seurat_obj, orig.ident == "DBCO")
exprMatrix <- GetAssayData(DBCO_obj, slot = "data")
cells_rankings <- AUCell_buildRankings(exprMatrix)
aucMaxRank<- 1000
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = aucMaxRank)
save(cells_AUC, file="DBCO_AUC.RData")

## 案例
library(Seurat)
library(AUCell)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
# 获取 AUCell 得分矩阵
seurat_obj <- readRDS("D:/File/document/单细胞测序代码/AUCell_tutorial/all_seurat_obj.rds") 
load("D:/File/document/单细胞测序代码/AUCell_tutorial/all_cells_AUC.RData")
# Step 1: Get AUCell matrix and convert to a data frame
# Step 1: Get AUCell matrix and convert to a data frame
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))

# Step 2: Add group info from Seurat object
auc_df$orig.ident <- seurat_obj$orig.ident[match(rownames(auc_df), colnames(seurat_obj))]
auc_df$celltype <- seurat_obj$celltype[match(rownames(auc_df), colnames(seurat_obj))]

# Step 3: Reshape the data into a long format
auc_long <- pivot_longer(auc_df, cols = -c(orig.ident, celltype),
                         names_to = "gene_set", values_to = "auc_score")

# Step 4: Calculate the mean AUCell for each combination of orig.ident, celltype, and gene_set
mean_df <- auc_long %>%
  group_by(orig.ident, celltype, gene_set) %>%
  summarise(mean_auc = mean(auc_score, na.rm = TRUE), .groups = "drop")

# Step 4: Calculate the sum of AUCell for each combination of orig.ident, celltype, and gene_set
sum_df <- auc_long %>%
  group_by(orig.ident, celltype, gene_set) %>%  # 按照 orig.ident, celltype 和 gene_set 进行分组
  summarise(sum_auc = sum(auc_score, na.rm = TRUE), .groups = "drop")  # 计算每组的 auc_score 总和，并去除缺失值

# Step 5: Get unique gene sets and sample identifiers
unique_gene_sets <- unique(mean_df$gene_set)
selected_gene_set_indices <- c(2, 3, 4, 6, 7, 9, 12)
selected_gene_sets <- unique_gene_sets[selected_gene_set_indices]
selected_celltypes <- c("T cells", "NK cells", "Monocytes", "Granulocytes", "Macrophages", "B cell")

# Step 6: Filter the data to select the specified gene sets and orig.ident
mean_df_selected <- mean_df %>% filter(gene_set %in% selected_gene_sets & 
                                        celltype %in% selected_celltypes)
sum_df_selected <- sum_df %>% filter(gene_set %in% selected_gene_sets & 
                                        celltype %in% selected_celltypes)
auc_long_selected <- auc_long %>% filter(gene_set %in% selected_gene_sets & 
                                          celltype %in% selected_celltypes)

# Step 7: Get the sample groups (orig.ident)
orig_idents <- unique(mean_df_selected$orig.ident)

# Step 8: Filter data for individual groups (e.g., DBCO, A0_EV, A1_EV, A2_EV)
dbco_data <- mean_df_selected %>% filter(orig.ident == "DBCO")
A0EV_data <- mean_df_selected %>% filter(orig.ident == "A0_EV")
A1EV_data <- mean_df_selected %>% filter(orig.ident == "A1_EV")
A2EV_data <- mean_df_selected %>% filter(orig.ident == "A2_EV")
stroke_data <- mean_df_selected %>% filter(orig.ident == "Stroke")
# Step 8: Filter data for individual groups (e.g., DBCO, A0_EV, A1_EV, A2_EV)
dbco_data <- sum_df_selected %>% filter(orig.ident == "DBCO")
A0EV_data <- sum_df_selected %>% filter(orig.ident == "A0_EV")
A1EV_data <- sum_df_selected %>% filter(orig.ident == "A1_EV")
A2EV_data <- sum_df_selected %>% filter(orig.ident == "A2_EV")
stroke_data <- sum_df_selected %>% filter(orig.ident == "Stroke")
# Step 9: Calculate the mean differences between experimental groups and DBCO (control group)
# Create a data frame for the mean difference
mean_diff_df <- A2EV_data %>%
  left_join(dbco_data, by = c("celltype", "gene_set"), suffix = c("_exp", "_dbco")) %>%
  mutate(mean_diff = mean_auc_exp - mean_auc_dbco)

sum_diff_df <- A2EV_data %>%
  left_join(dbco_data, by = c("celltype", "gene_set"), suffix = c("_exp", "_dbco")) %>%
  mutate(sum_diff = sum_auc_exp - sum_auc_dbco)
# Step 10: Reshape the data into long format if not already done
# 确保 auc_long_selected 是长格式数据
selected_groups <- c("DBCO", "A2_EV")  # 只保留对照组和实验组
auc_long_selected_filtered <- auc_long_selected%>%filter(orig.ident %in% selected_groups)
stat_df_cells <- auc_long_selected_filtered %>%
  group_by(celltype, gene_set) %>%
  summarise(
    p_value = kruskal.test(auc_score ~ orig.ident)$p.value,  # 对每个 celltype 和 gene_set 组合进行 Kruskal-Wallis 检验
    .groups = "drop"
  ) %>%
  mutate(log10_p = ifelse(p_value == 0, 100, -log10(p_value)))  # 将 p 值转换为 -log10(p)，以便更直观地显示
stat_df_cells

# Step 11: Merge the mean differences and statistical significance data
plot_df <- left_join(mean_diff_df, stat_df_cells, by = c("celltype", "gene_set"))
plot_df <- left_join(sum_diff_df, stat_df_cells, by = c("celltype", "gene_set"))
# 将 log10_p 中的 Inf 替换为 300，其他值保持不变
plot_df <- plot_df %>%
  mutate(log10_p = ifelse(is.infinite(log10_p), 300, log10_p))

# Create the plot
p <- ggplot(plot_df, aes(x = gene_set, y = celltype, color = mean_diff, size =log10_p)) +
  geom_point(stroke = 0) +
  scale_color_gradient2(low = "#0073C2",mid ="#ffcddd",high = "#ff1a53", midpoint = 0) +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Functional gene set activity difference for selected groups vs DBCO",
    x = "Gene Set",
    y = "Cell Type",
    color = "Difference in Mean AUCell",
    size = "Abs Difference in -log10(p)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

# Display the plot
print(p)
# 保存图像，大小为 10x10 英寸
ggsave("kruskal_test_results.png", plot = p, width = 10, height = 10, units = "in", dpi = 300)

library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)  # 用于 pivot_longer

# 读取数据
seurat_obj <- readRDS("D:/File/document/单细胞测序代码/AUCell_tutorial/all_seurat_obj.rds") 
load("D:/File/document/单细胞测序代码/AUCell_tutorial/all_cells_AUC.RData")  # 加载 cells_AUC

# 获取 AUCell 得分矩阵并转为数据框
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))

# 添加元信息（group 和 celltype）
auc_df$orig.ident <- seurat_obj$orig.ident[match(rownames(auc_df), colnames(seurat_obj))]
auc_df$celltype <- seurat_obj$celltype[match(rownames(auc_df), colnames(seurat_obj))]

# 转为长表形式
auc_long <- pivot_longer(auc_df, cols = -c(orig.ident, celltype),
                         names_to = "gene_set", values_to = "auc_score")

# 添加 Group 列：将 orig.ident 中包含 "EV" 的标记为 EV，其他为 DBCO
auc_long$Group <- ifelse(grepl("EV", auc_long$orig.ident), "EV", "DBCO")


# 过滤出需要的数据
plot_df <- auc_long %>%
  filter(
    gene_set %in% c("myelination", "positive_regulation_of_angiogenesis"),
    celltype %in% c("T cells", "NK cells", "Monocytes", "Granulocytes", "Macrophages", "B cell"),
    orig.ident %in% c("A0_EV", "A1_EV", "A2_EV", "DBCO")
  )

# 设置颜色方案（为 A0_EV, A1_EV, A2_EV 和 DBCO 设置不同的颜色）
color_map <- c("A0_EV" = "#FF6347", "A1_EV" = "#4682B4", "A2_EV" = "#32CD32", "DBCO" = "#FFD700")

# 分别绘图并保存
gene_sets_to_plot <- c("myelination", "positive_regulation_of_angiogenesis")

for (gs in gene_sets_to_plot) {
  p <- ggplot(plot_df %>% filter(gene_set == gs),
              aes(x = auc_score, y = celltype, fill = orig.ident, color = orig.ident)) +
    geom_density_ridges(alpha = 0.4, scale = 1.0, rel_min_height = 0.01, size = 0.6) +  # 透明度调整为0.4，外围线条大小为0.6
    theme_minimal(base_size = 13) +
    # 使用定义的颜色映射
    scale_fill_manual(values = color_map) +
    scale_color_manual(values = color_map) +
    labs(title = paste0("Gene Set: ", gs),
         x = "AUCell Score", y = "Cell Type") +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none"
    )
  
  print(p)  # 显示图
  
  # 保存为 PNG 文件
  ggsave(filename = paste0("AUCell_", gs, ".png"), plot = p, width = 7, height = 5, dpi = 300)
}
