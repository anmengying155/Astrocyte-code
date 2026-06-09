# ============================================================================
# Macrophage Subgroup Analysis
# Purpose: Analyze MF1/MF2/MF3 macrophage subgroups with marker genes
# Input: RDS file with macrophage/neutrophil data
# Output: Violin plots, bar charts, statistics
# ============================================================================

library(Seurat)
library(ggplot2)

# Define paths relative to project root
SCRIPT_DIR <- dirname(parent.frame(2)$ofile)
PROJECT_ROOT <- ifelse(is.null(SCRIPT_DIR), ".", dirname(SCRIPT_DIR))

DATA_DIR <- file.path(PROJECT_ROOT, "..", "..", "Seq2_3_Comparison", "summary (2)", "summary", "2_Expression_result")
RDS_FILE <- file.path(DATA_DIR, "A2_EV_macrophage_neutrophil.rds")
OUTPUT_DIR <- DATA_DIR

seurat_obj <- readRDS(RDS_FILE)
cat("Cells:", ncol(seurat_obj), "\n")
cat("Original Idents:", levels(seurat_obj@active.ident), "\n\n")

cluster_mapping <- c(
  "0" = "MF1",
  "1" = "MF2",
  "2" = "MF3"
)

current_ids <- as.character(Idents(seurat_obj))
new_ids <- cluster_mapping[current_ids]
Idents(seurat_obj) <- factor(new_ids, levels = c("MF1", "MF2", "MF3"))

cat("New Idents:", levels(seurat_obj@active.ident), "\n\n")

markers <- list(
  MF1 = c("Col6a3", "Hsd11b2", "Ogn", "Col6a4", "Igfbp5"),
  MF2 = c("H2-Aa", "H2-Eb1", "H2-Ab1", "C3", "Saa3"),
  MF3 = c("Igkc", "Igha", "Igfbp7", "Dpt", "Arap1")
)

cat("=== Violin plots for marker genes ===\n")
for (mf_name in names(markers)) {
  genes <- markers[[mf_name]]
  genes_in_data <- genes[genes %in% rownames(seurat_obj)]
  if (length(genes_in_data) == 0) {
    cat(mf_name, "- no marker genes found\n")
    next
  }

  plot_data_list <- list()
  for (gene in genes_in_data) {
    expr_data <- as.numeric(seurat_obj@assays$RNA@data[gene, ])
    plot_data_list[[gene]] <- data.frame(
      Gene = gene,
      Expression = expr_data,
      CellType = Idents(seurat_obj),
      Sample = seurat_obj$Sample
    )
  }
  plot_df <- do.call(rbind, plot_data_list)

  p <- ggplot(plot_df, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_violin(trim = FALSE) +
    scale_fill_manual(values = c("MF1" = "#1f77b4", "MF2" = "#ff7f0e", "MF3" = "#2ca02c")) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(mf_name, "Marker Genes"), x = "Cell Type", y = "Expression") +
    facet_wrap(~Gene, ncol = 5, scales = "free_y")

  pngfilename <- file.path(OUTPUT_DIR, paste0(mf_name, "_markers_violin.png"))
  png(pngfilename, width = 15, height = 4, units = "in", res = 300)
  print(p)
  dev.off()
  cat("Saved:", basename(pngfilename), "\n")
}

cat("\n=== Sample statistics ===\n")
stat_df <- data.frame(
  Sample = seurat_obj$Sample,
  CellType = Idents(seurat_obj)
)
stat_table <- table(stat_df$Sample, stat_df$CellType)
stat_percent <- prop.table(stat_table, margin = 1) * 100

cat("\nCell counts:\n")
print(stat_table)
cat("\nCell percentages:\n")
print(round(stat_percent, 2))

stat_summary <- as.data.frame(stat_percent)
colnames(stat_summary) <- c("Sample", "CellType", "Percentage")
stat_summary$Count <- as.vector(stat_table)

write.csv(stat_summary, file.path(OUTPUT_DIR, "MF_sample_statistics.csv"), row.names = FALSE)
cat("\nSaved: MF_sample_statistics.csv\n")

p <- ggplot(stat_summary, aes(x = Sample, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("MF1" = "#1f77b4", "MF2" = "#ff7f0e", "MF3" = "#2ca02c")) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), position = position_dodge(0.9), vjust = -0.5, size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Macrophage Subgroup Distribution by Sample", x = "Sample", y = "Percentage (%)", fill = "Cell Type")

pngfilename <- file.path(OUTPUT_DIR, "MF_sample_distribution.png")
png(pngfilename, width = 8, height = 6, units = "in", res = 300)
print(p)
dev.off()
cat("Saved:", basename(pngfilename), "\n")

cat("\nDone!\n")
