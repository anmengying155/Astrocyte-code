# -------------------------------
# 加载必要的库
# -------------------------------
library(SingleR)
library(celldex)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(openxlsx)

load("D:/File/document/单细胞测序代码/singleR注释库/sub_巨噬+单核_cells_去批次_bimod_Filter_Basicinfo.rData")
# 查看基本信息
head(pbmc_filt)
celltype_info <- read.xlsx("D:/File/document/单细胞测序代码/cell_annot (6-Macro).xlsx")
head(celltype_info)
seurat_obj <- pbmc_filt
seurat_obj$celltype <- celltype_info$cells[match(seurat_obj$seurat_clusters, celltype_info$cluster)]
head(seurat_obj)
DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", pt.size = 3) + 
  ggtitle("seurat_umap") + coord_fixed()

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
aucMaxRank<- 1436
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = aucMaxRank)
save(cells_AUC, file = "巨噬再分群_AUC.RData")



## 案例
library(Seurat)
library(AUCell)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggridges)
# Step 1: Get AUCell matrix and convert to a data frame
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))
auc_df$ celltype <- seurat_obj$celltype[match(rownames(auc_df), colnames(seurat_obj))]
auc_df$orig.ident <- seurat_obj$orig.ident[match(rownames(auc_df), colnames(seurat_obj))]

# Step 3: Reshape the data into a long format
auc_long <- pivot_longer(auc_df, cols = -c(orig.ident,celltype),
                         names_to = "gene_set", values_to = "auc_score")

# Step 4: Calculate the mean AUCell for each combination of celltype, and gene_set
mean_df <- auc_long %>%
  group_by(celltype, gene_set,orig.ident) %>%
  summarise(mean_auc = mean(auc_score, na.rm = TRUE), .groups = "drop")

# Step 5: Get unique gene sets and sample identifiers
unique_gene_sets <- unique(auc_long$gene_set)
selected_gene_set_indices <- c(6)
selected_gene_sets <- unique_gene_sets[selected_gene_set_indices]
unique_celltypes <- unique(auc_long$celltype)
selected_celltypes_indices <- c(1,4,6,7,8,10)
selected_celltypes<- unique_celltypes[selected_celltypes_indices]
unique_orig.ident <- unique(auc_long$orig.ident)
selected_orig.ident_indices <- c(1,2,3)
selected_orig.ident<- unique_orig.ident[selected_orig.ident_indices]
# Step 6: Filter the data to select the specified gene sets and orig.ident
auc_long_selected <- auc_long %>% filter(gene_set %in% selected_gene_sets & celltype %in% selected_celltypes & orig.ident %in% selected_orig.ident )
mean_df_selected <- mean_df%>% filter(gene_set %in% selected_gene_sets & celltype %in% selected_celltypes & orig.ident %in% selected_orig.ident )
# Create a data frame for the mean difference
write.csv(as.data.frame(mean_df_selected), "mean_df_selected.csv")

# 根据中位数排序
label_order <- auc_long_selected %>%
  group_by(celltype) %>%
  summarize(median_auc = median(auc_score)) %>%
  arrange(desc(median_auc)) %>%
  pull(labels)

# 应用排序
auc_long_selected$labels <- factor(auc_long_selected$labels, levels = label_order)

color_map <- c("A0_EV" = "#FF6347", "A1_EV" = "#4682B4", "A2_EV" = "#32CD32")

    p <- ggplot(auc_long_selected,
              aes(x = auc_score, y = celltype , fill = orig.ident, color = orig.ident)) +
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




  # 2. 数据转换为宽格式
heatmap_data <- mean_df_selected %>%
  pivot_wider(names_from = celltype, values_from = mean_auc)

# 3. 将 gene_set 设置为行名，并转换为 data.frame
heatmap_data <- as.data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data$orig.ident
heatmap_data <- heatmap_data %>% dplyr::select(-orig.ident)

# 4. 标准化数据（行或列）#00ff26#cbf7d1#cbf7e8##1a1a53ea#c4ffe6#bceed9
heatmap_data <- t(scale(t(heatmap_data)))  # 按行标准化（每个基因）

png("巨噬分群heatmap_output.png", width = 1000, height = 800)  # 设置文件名和图形的宽高

# 绘制热图，使用从 royalblue 到 firebrick3 的颜色渐变
pheatmap(heatmap_data, 
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
write.csv(heatmap_data, "巨噬分群heatmap_data.csv", row.names = TRUE)
  