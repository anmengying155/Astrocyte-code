library(Seurat)
library(AUCell)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
# 获取 AUCell 得分矩阵
load("D:/File/document/单细胞测序代码/AUCell_tutorial/all_cells_AUC.RData")
auc_mat <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_mat))
unique(seurat_obj$orig.ident)
# 添加 group info
auc_df$orig.ident <- seurat_obj$orig.ident[match(rownames(auc_df), colnames(seurat_obj))]
auc_df$celltype <- seurat_obj$celltype[match(rownames(auc_df), colnames(seurat_obj))]

# 转为长表
auc_long <- pivot_longer(auc_df, cols = -c(orig.ident, celltype),
                         names_to = "gene_set", values_to = "auc_score")

# 计算每个组合的 AUCell 均值
mean_df <- auc_long %>%
  group_by(orig.ident, celltype, gene_set) %>%
  summarise(mean_auc = mean(auc_score, na.rm = TRUE), .groups = "drop")

# 计算每种 celltype + gene_set 下，orig.ident 组之间是否显著不同
# 使用 Kruskal-Wallis 检验（适合多组）
stat_df <- auc_long %>%
  group_by(celltype, gene_set) %>%
  summarise(
    p_value = kruskal.test(auc_score ~ orig.ident)$p.value,
    .groups = "drop"
  ) %>%
  mutate(log10_p = -log10(p_value))
# 合并均值和显著性信息
plot_df <- left_join(mean_df, stat_df, by = c("celltype", "gene_set"))
save(plot_df, file = "plot_df.RData")

# 获取 gene_set 中的唯一值和 orig.ident（样本或组别）
unique_gene_sets <- unique(plot_df$gene_set)
unique_orig_idents <- unique(plot_df$orig.ident)
# 假设你要选择特定的基因集合
selected_gene_sets <- unique_gene_sets[selected_gene_set_indices]
plot_df_selected <- plot_df %>% filter(gene_set %in% selected_gene_sets)

# 获取所有的 orig.ident（样本或组别）
orig_idents <- unique(plot_df_selected$orig.ident)

# 排除 DBCO 组
dbco_data <- plot_df_selected %>% filter(orig.ident == "DBCO")

# 循环保存每个 orig.ident 的单独图像，减去 DBCO 组的结果
for (ident in orig_idents) {
  # 过滤数据，只保留当前 orig.ident 的子集
  plot_data <- plot_df_selected %>% filter(orig.ident == ident)
  # 计算每个组与 DBCO 组的差值
  plot_data <- plot_data %>%
    left_join(dbco_data %>% dplyr::select(gene_set, celltype, mean_auc, log10_p), by = c("gene_set", "celltype"), suffix = c("", "_dbco"))
    mutate(
      diff_mean_auc = mean_auc - mean_auc_dbco,
      diff_log10_p = log10_p - log10_p_dbco
    ) %>%
    # 只选择差值大于零的数据
    filter(diff_mean_auc > 0, diff_log10_p > 0)
  
  # 如果差值大于零的数据存在，绘制 dotplot 图
  if (nrow(plot_data) > 0) {
    p <- ggplot(plot_data, aes(x = gene_set, y = celltype, color = diff_mean_auc, size = diff_log10_p)) +
      geom_point() +
      scale_color_gradient(low = "#ffebeb", high = "#ff6767") +
      theme_minimal(base_size = 13) +
      labs(
        title = paste("Functional gene set activity difference for", ident, "vs DBCO"),
        x = "Gene Set",
        y = "Cell Type",
        color = "Difference in Mean AUCell",
        size = "Difference in -log10(p)"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12)
      )
    
    # 保存为 SVG 文件
    ggsave(
      filename = paste0("dotplot_", ident, "_vs_DBCO_selected_diff_positive.svg"),
      plot = p,
      device = "svg",
      width = 10,
      height = 6
    )
  }
}
