# -------------------------------
# 3. 读取10X数据
# -------------------------------
data_path <- "F:/测序结果/第三次-星胶EV测序/Summary (1)/Summary/2_Expression_result/filtered_feature_bc_matrix" 
source_data <- Read10X(data.dir = data_path)
seurat_obj <- CreateSeuratObject(counts = source_data, project = "10X_Data", min.cells = 3, min.features = 200)
# 对 Seurat 对象中的数据进行标准化
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
head(seurat_obj)
head(orig_ident_info)
cluster_info <- read.csv("f:/测序结果/第三次-星胶EV测序/Summary (1)/Summary/2_Expression_result/Loupe_Browser_File/seurat_loupe_cell_umap_cluster.csv")
orig_ident_info <- read.csv("f:/测序结果/第三次-星胶EV测序/Summary (1)/Summary/2_Expression_result/orig.ident.csv")
barcode_info$Barcode <- gsub(" ", "", cluster_info$Cell)
barcode_info$seurat_clusters <- gsub(" ", "", cluster_info$Cluster)
barcode_info$orig.ident <- gsub(" ", "", orig_ident_info$orig.ident)

seurat_obj$orig.ident <- barcode_info$orig.ident[match(Cells(seurat_obj), barcode_info$Barcode)]
seurat_obj$seurat_clusters <- barcode_info$seurat_clusters[match(Cells(seurat_obj), barcode_info$Barcode)]

# 取出非 NA 的细胞名
valid_cells <- colnames(seurat_obj)[!is.na(seurat_obj$orig.ident)]
# 用细胞名子集化 Seurat 对象
seurat_obj <- subset(seurat_obj, cells = valid_cells)
# 查看 Seurat 对象的所有元数据列
colnames(seurat_obj@meta.data)
# 4. 查看添加后的 Seurat 对象的元数据
head(seurat_obj@meta.data)
table(is.na(seurat_obj$orig.ident))  # 看看多少 NA
unique(seurat_obj$orig.ident)
seurat_obj <- NormalizeData(seurat_obj) 
# 保存 Seurat 对象
saveRDS(seurat_obj, file = "A0-A2_seurat_obj.rds")
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
#write.csv(cell_counts_MR, file = "cell_counts_MR.csv", row.names = TRUE)

endothelial_cells <- subset(seurat_obj, subset = MR_labels == "Endothelial cells")
# 绘制内皮细胞的 UMAP 图
dimplot_endothelial_cells<- DimPlot(endothelial_cells, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of Endothelial Cells")

# 创建子集，选择 seurat_clusters 为 "1", "21", "23" 的细胞
sub_seurat_obj <- subset(seurat_obj, subset = seurat_clusters %in% c("1", "21", "23","26"))
# 绘制内皮细胞的 UMAP 图
dimplot_sub_seurat_obj <- DimPlot(sub_seurat_obj, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of Endothelial Cells")


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
exprMatrix <- GetAssayData(sub_seurat_obj, slot = "data")
cells_rankings <- AUCell_buildRankings(exprMatrix)
aucMaxRank<- 1033
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = aucMaxRank)
save(cells_AUC, file="A0-A2_endothelial_cells_AUC.RData")

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
# Step 1: Get AUCell matrix and convert to a data frame
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))

# Step 2: Add group info from Seurat object
auc_df$orig.ident <- seurat_obj$orig.ident[match(rownames(auc_df), colnames(seurat_obj))]
# Step 3: Reshape the data into a long format
auc_long <- pivot_longer(auc_df, cols = -c(orig.ident),
                         names_to = "gene_set", values_to = "auc_score")

# Step 5: Get unique gene sets and sample identifiers
unique_gene_sets <- unique(auc_long$gene_set)
unique_orig.ident <- unique(auc_long$orig.ident)
selected_gene_set_indices <- c(6)
selected_gene_sets <- unique_gene_sets[selected_gene_set_indices]
# Step 6: Filter the data to select the specified gene sets and orig.ident
auc_long_selected <- auc_long %>% filter(gene_set %in% selected_gene_sets )



color_map <- c("MO" = "#FF6347", "REP" = "#4682B4", "NT" = "#32CD32")

    p <- ggplot(auc_long_selected,
              aes(x = auc_score, y = orig.ident , fill = orig.ident, color = orig.ident)) +
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
  
