library(ggplot2)
library(dplyr)
library(readr)

miRNA_data <- read_csv('d:/单细胞测序/Original_Code/05_ProcessedData/Expression_Data/miRNA_list.csv')

miRNA_data <- miRNA_data %>%
  mutate(log2FC = log2(`fold_change(A2_EV/A1_EV)`)) %>%
  mutate(log10_p = -log10(`pvalue(chi_square_2x2)`)) %>%
  mutate(log2FC = ifelse(is.infinite(log2FC) & log2FC > 0, 10, log2FC)) %>%
  mutate(log2FC = ifelse(is.infinite(log2FC) & log2FC < 0, -10, log2FC)) %>%
  mutate(log10_p = ifelse(is.infinite(log10_p), 80, log10_p)) %>%
  mutate(log10_p = ifelse(is.na(log10_p), 0, log10_p)) %>%
  mutate(
    significance = case_when(
      abs(log2FC) >= 1 & log10_p >= 1.3 ~ '显著',
      TRUE ~ '不显著'
    )
  ) %>%
  mutate(
    direction = case_when(
      log2FC >= 1 & log10_p >= 1.3 ~ '上调',
      log2FC <= -1 & log10_p >= 1.3 ~ '下调',
      TRUE ~ '无显著变化'
    )
  )

sig_count <- miRNA_data %>% filter(significance == '显著') %>% nrow()
up_count <- miRNA_data %>% filter(direction == '上调') %>% nrow()
down_count <- miRNA_data %>% filter(direction == '下调') %>% nrow()

cat(sprintf('总 miRNA 数: %d\n', nrow(miRNA_data)))
cat(sprintf('显著差异 miRNA 数: %d\n', sig_count))
cat(sprintf('上调 miRNA 数: %d\n', up_count))
cat(sprintf('下调 miRNA 数: %d\n', down_count))

top_up <- miRNA_data %>% filter(direction == '上调') %>% arrange(desc(log2FC)) %>% head(5)
top_down <- miRNA_data %>% filter(direction == '下调') %>% arrange(log2FC) %>% head(5)

cat('\nTop 5 上调 miRNA:\n')
print(top_up %>% select(miR_name, `fold_change(A2_EV/A1_EV)`, `pvalue(chi_square_2x2)`, log2FC))

cat('\nTop 5 下调 miRNA:\n')
print(top_down %>% select(miR_name, `fold_change(A2_EV/A1_EV)`, `pvalue(chi_square_2x2)`, log2FC))

plot_df <- miRNA_data %>% filter(abs(log2FC) <= 12 & log10_p <= 85)

p <- ggplot(plot_df, aes(x = log2FC, y = log10_p, color = direction)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'gray50', size = 0.5) +
  geom_hline(yintercept = 1.3, linetype = 'dashed', color = 'gray50', size = 0.5) +
  scale_color_manual(values = c(
    '上调' = '#E41A1C',
    '下调' = '#377EB8',
    '无显著变化' = '#999999'
  )) +
  labs(
    x = 'log2(Fold Change)',
    y = '-log10(p-value)',
    title = 'miRNA 差异表达火山图',
    subtitle = sprintf('A2_EV vs A1_EV | 显著差异: %d (上调: %d, 下调: %d)', sig_count, up_count, down_count),
    color = '差异表达'
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = 'gray90'),
    panel.grid.minor = element_line(color = 'gray95')
  )

pdf('d:/单细胞测序/scRNAseq_Code/Main_Scripts/Output/miRNA_volcano_plot.pdf', width = 10, height = 8)
print(p)
dev.off()

png('d:/单细胞测序/scRNAseq_Code/Main_Scripts/Output/miRNA_volcano_plot.png', width = 1000, height = 800, res = 100)
print(p)
dev.off()

cat('\n火山图已保存至:\n')
cat('  - miRNA_volcano_plot.pdf\n')
cat('  - miRNA_volcano_plot.png\n')
cat('输出目录: d:/单细胞测序/scRNAseq_Code/Main_Scripts/Output/\n')