#案例
library(Seurat)
library(AUCell)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
# 获取 AUCell 得分矩阵
seurat_obj <- readRDS("D:/File/document/单细胞测序代码/AUCell_tutorial/all_seurat_obj.rds") 
load("D:/File/document/单细胞测序代码/A2_EV_AUC.RData")
# Step 1: Get AUCell matrix and convert to a data frame
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))
# Step 2: Add group info from Seurat object
auc_df$orig.ident <- A2_EV_obj$orig.ident[match(rownames(auc_df), colnames(A2_EV_obj))]
auc_df$celltype <- A2_EV_obj$celltype[match(rownames(auc_df), colnames(A2_EV_obj))]

# Step 3: Reshape the data into a long format
auc_long1 <- pivot_longer(auc_df, cols = -c(orig.ident, celltype),
                         names_to = "gene_set", values_to = "auc_score")

# Step 4: Calculate the mean AUCell for each combination of orig.ident, celltype, and gene_set
mean_df1 <- auc_long1 %>%
  group_by(orig.ident, celltype, gene_set) %>%
  summarise(mean_auc = mean(auc_score, na.rm = TRUE), .groups = "drop")


load("D:/File/document/单细胞测序代码/DBCO_AUC.RData")
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))
# Step 2: Add group info from Seurat object
auc_df$orig.ident <- DBCO_obj$orig.ident[match(rownames(auc_df), colnames(DBCO_obj))]
auc_df$celltype <- DBCO_obj$celltype[match(rownames(auc_df), colnames(DBCO_obj))]
# Step 3: Reshape the data into a long format
auc_long2 <- pivot_longer(auc_df, cols = -c(orig.ident, celltype),
                         names_to = "gene_set", values_to = "auc_score")

# Step 4: Calculate the mean AUCell for each combination of orig.ident, celltype, and gene_set
mean_df2 <- auc_long2 %>%
  group_by(orig.ident, celltype, gene_set) %>%
  summarise(mean_auc = mean(auc_score, na.rm = TRUE), .groups = "drop")

# 合并 mean_df1 和 mean_df2
mean_df <- left_join(mean_df1, mean_df2, by = c("celltype", "gene_set"), suffix = c("_A2_EV", "_DBCO"))
# 计算 A2_EV 和 DBCO 的 AUCell 值差异
mean_df <- mean_df %>%
  mutate(diff_mean_auc = mean_auc_A2_EV - mean_auc_DBCO)
# 查看结果
head(mean_df)
# Step 5: Get unique gene sets and sample identifiers
unique_gene_sets <- unique(mean_df$gene_set)
selected_gene_set_indices <- c(2, 3, 4, 6, 7, 9, 12)
selected_gene_sets <- unique_gene_sets[selected_gene_set_indices]
#selected_celltypes <- c("T cells", "NK cells", "Monocytes", "Granulocytes", "Macrophages", "B cell")

# Step 6: Filter the data to select the specified gene sets and orig.ident
mean_df_selected <- mean_df %>% filter(gene_set %in% selected_gene_sets)
auc_long1_selected <- auc_long1 %>% filter(gene_set %in% selected_gene_sets )
auc_long2_selected <- auc_long2 %>% filter(gene_set %in% selected_gene_sets )

# Step 1: 添加 group 列（标识不同数据集）
auc_long1_selected$group <- "Group1"  # 添加 Group1 标签
auc_long2_selected$group <- "Group2"  # 添加 Group2 标签

# Step 2: 合并两个数据集
auc_combined <- bind_rows(auc_long1_selected, auc_long2_selected)
# Step 3: 进行 Kruskal-Wallis 检验
stat_df_cells <- auc_combined %>%
  group_by(celltype, gene_set) %>%
  summarise(
    p_value = kruskal.test(auc_score ~ group)$p.value,  # 对每个 celltype 和 gene_set 组合进行 Kruskal-Wallis 检验
    .groups = "drop"
  ) %>%
  mutate(log10_p = ifelse(p_value == 0, 100, -log10(p_value)))  # 将 p 值转换为 -log10(p)，以便更直观地显示

# 查看检验结果
stat_df_cells

# Step 11: Merge the mean differences and statistical significance data
plot_df <- left_join(mean_df_selected, stat_df_cells, by = c("celltype", "gene_set"))
# 将 log10_p 中的 Inf 替换为 300，其他值保持不变
plot_df <- plot_df %>%
  mutate(log10_p = ifelse(is.infinite(log10_p), 300, log10_p))

# 筛选出 diff_mean_auc >= 0 的数据
plot_df_filtered <- plot_df %>% filter(diff_mean_auc >=-0.03)

# 创建图表
p <- ggplot(plot_df_filtered, aes(x = gene_set, y = celltype, color = diff_mean_auc, size = log10_p)) +
  geom_point(stroke = 0) +
  scale_color_gradient2(low = "#0073C2", mid = "#ffcddd", high = "#ff1a53", midpoint = 0) +
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
# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)


# 1. 去除 orig.ident 列
# Step 5: Get unique gene sets and sample identifiers
unique_gene_sets <- unique(mean_df$gene_set)
selected_gene_set_indices <- c(2, 3, 4, 6, 7, 9, 12)
selected_gene_sets <- unique_gene_sets[selected_gene_set_indices]
#selected_celltypes <- c("T cells", "NK cells", "Monocytes", "Granulocytes", "Macrophages", "B cell")

# Step 6: Filter the data to select the specified gene sets and orig.ident
mean_df_selected <- mean_df1 %>% filter(gene_set %in% selected_gene_sets)
mean_df_clean <- mean_df_selected %>%
  dplyr::select(-orig.ident)

# 2. 数据转换为宽格式
heatmap_data <- mean_df_clean %>%
  pivot_wider(names_from = celltype, values_from = mean_auc)

# 3. 将 gene_set 设置为行名，并转换为 data.frame
heatmap_data <- as.data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data$gene_set
heatmap_data <- heatmap_data %>% dplyr::select(-gene_set)

# 4. 标准化数据（行或列）#00ff26#cbf7d1#cbf7e8##1a1a53ea#c4ffe6#bceed9
heatmap_data <- t(scale(t(heatmap_data)))  # 按行标准化（每个基因）

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
write.csv(heatmap_data, "heatmap_data.csv", row.names = TRUE)