# ============================================================================
# Selected miRNA Target Analysis
# Purpose: Analyze specific miRNA targets with higher expression threshold
# Input: RDS file with macrophage/neutrophil data, miRNA full dataset
# Output: Heatmap and violin plots for selected targets
# ============================================================================

library(Seurat)
library(readxl)
library(ggplot2)
library(pheatmap)

# Define paths relative to project root
SCRIPT_DIR <- dirname(parent.frame(2)$ofile)
PROJECT_ROOT <- ifelse(is.null(SCRIPT_DIR), ".", dirname(SCRIPT_DIR))

DATA_DIR <- file.path(PROJECT_ROOT, "..", "..", "Seq2_3_Comparison", "summary (2)", "summary", "2_Expression_result")
RDS_FILE <- file.path(DATA_DIR, "A2_EV_macrophage_neutrophil.rds")
MICRORNA_FILE <- file.path(PROJECT_ROOT, "..", "Data", "miRNA_full_dataset_source.xlsx")
EXPR_THRESHOLD <- 0.5

seurat_obj <- readRDS(RDS_FILE)
microrna_df <- read_excel(MICRORNA_FILE)

cat("Cells in dataset:", ncol(seurat_obj), "\n\n")

target_mirnas <- c("rno-miR-28-3p_R+1", "rno-miR-93-5p", "rno-miR-351-5p_R-2", "rno-miR-450a-5p", "rno-miR-542-3p")
cell_types <- c("Macrophage", "Neutrophil")

all_targets <- c()
for (mirna in target_mirnas) {
  targets <- microrna_df$Symbol[microrna_df$`miRNA ID` == mirna]
  all_targets <- unique(c(all_targets, targets))
}
cat("Total unique targets from selected miRNAs:", length(all_targets), "\n")

all_targets <- all_targets[all_targets %in% rownames(seurat_obj)]
cat("Targets found in scRNA-seq data:", length(all_targets), "\n")

plot_data_matrix <- matrix(0, nrow = length(all_targets), ncol = length(cell_types))
rownames(plot_data_matrix) <- all_targets
colnames(plot_data_matrix) <- cell_types

for (cell_type in cell_types) {
  seurat_subset <- seurat_obj[, seurat_obj$CellType == cell_type]
  layer_name <- Layers(seurat_subset)[1]
  expr_matrix <- LayerData(seurat_subset, layer = layer_name)[all_targets, , drop = FALSE]
  expr_means <- rowMeans(expr_matrix)
  plot_data_matrix[, cell_type] <- expr_means
}

max_expr <- apply(plot_data_matrix, 1, max)
filtered_targets <- names(max_expr[max_expr >= EXPR_THRESHOLD])
cat("Targets with expression >= 0.5:", length(filtered_targets), "\n\n")

plot_data_matrix <- plot_data_matrix[filtered_targets, , drop = FALSE]
plot_data_matrix <- plot_data_matrix[order(rowMeans(plot_data_matrix), decreasing = TRUE), , drop = FALSE]

max_expr_all <- max(plot_data_matrix)
break_vals <- seq(0, ceiling(max_expr_all * 10) / 10, length.out = 20)
color_scheme <- colorRampPalette(c("blue", "red"))(length(break_vals) - 1)

cat("=== Generating combined heatmap ===\n")
pngfilename <- file.path(DATA_DIR, "A2_EV_selected_targets_heatmap.png")
png(pngfilename, width = 6, height = max(6, nrow(plot_data_matrix) * 0.3), units = "in", res = 300)
pheatmap(plot_data_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         legend = TRUE,
         breaks = break_vals,
         color = color_scheme,
         main = "Selected Targets Expression (Mean)")
dev.off()
cat("Saved:", basename(pngfilename), "\n\n")

cat("=== Generating violin plots ===\n")
for (target in filtered_targets) {
  plot_data_list_ct <- list()

  for (cell_type in cell_types) {
    seurat_subset <- seurat_obj[, seurat_obj$CellType == cell_type]
    layer_name <- Layers(seurat_subset)[1]
    expr <- LayerData(seurat_subset, layer = layer_name)[target, ]
    expr_df <- data.frame(Expression = as.numeric(expr), CellType = cell_type)
    plot_data_list_ct[[cell_type]] <- expr_df
  }

  combined_data <- do.call(rbind, plot_data_list_ct)

  p <- ggplot(combined_data, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_violin(trim = FALSE) +
    geom_jitter(size = 0.1, alpha = 0.3) +
    scale_fill_manual(values = c("Macrophage" = "blue", "Neutrophil" = "red")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = target, x = "Cell Type", y = "Expression")

  pngfilename <- file.path(DATA_DIR, paste0("A2_EV_", target, "_violin.png"))
  ggsave(pngfilename, plot = p, width = 6, height = 5, dpi = 300)
  cat("Saved:", basename(pngfilename), "\n")
}

cat("\nDone!\n")
