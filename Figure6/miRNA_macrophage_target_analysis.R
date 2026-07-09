# ============================================================================
# Macrophage miRNA Target Analysis
# Purpose: Analyze miRNA targets in Macrophage vs Neutrophil cells
# Input: RDS files with cell data, miRNA target predictions
# Output: Violin plots, heatmaps, bubble plots, CSV data
# ============================================================================

library(Seurat)
library(readxl)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)
library(svglite)

# Define paths relative to project root
SCRIPT_DIR <- dirname(parent.frame(2)$ofile)
PROJECT_ROOT <- ifelse(is.null(SCRIPT_DIR), ".", dirname(SCRIPT_DIR))

# ===================== Configuration =====================
DATA_DIR <- file.path(PROJECT_ROOT, "..", "..", "Seq2_3_Comparison", "summary (2)", "summary", "2_Expression_result")
RDS_FILES <- c(
  file.path(DATA_DIR, "A0_EV_macrophage_neutrophil.rds"),
  file.path(DATA_DIR, "A1_EV_macrophage_neutrophil.rds"),
  file.path(DATA_DIR, "A2_EV_macrophage_neutrophil.rds")
)
MICRORNA_FILE <- file.path(PROJECT_ROOT, "..", "Data", "miRNA_target_predictions.csv")
EXPR_THRESHOLD <- 0

microrna_df <- read.csv(MICRORNA_FILE)
target_mirnas <- unique(microrna_df$miRNA.ID)
cell_types <- c("Macrophage", "Neutrophil")

for (RDS_FILE in RDS_FILES) {
  cat("\n========== Processing:", basename(RDS_FILE), "==========\n")
  seurat_obj <- readRDS(RDS_FILE)
  cat("Cells in dataset:", ncol(seurat_obj), "\n")
  cat("Target miRNAs:", length(target_mirnas), "\n\n")

  all_csv_data <- data.frame(stringsAsFactors = FALSE)
  max_expr_all <- 0
  plot_data_list <- list()
  
  sample_name <- gsub("_macrophage_neutrophil\\.rds$", "", basename(RDS_FILE))

# ===================== 遍历miRNA计算平均表达量 =====================
for (mirna in target_mirnas) {
  targets <- microrna_df$Symbol[microrna_df$miRNA.ID == mirna]
  targets <- targets[targets %in% rownames(seurat_obj)]
  if (length(targets) == 0) next
  
  means_list <- c()
  for (cell_type in cell_types) {
    seurat_subset <- seurat_obj[, seurat_obj$CellType == cell_type]
    layer_name <- Layers(seurat_subset)[1]
    plot_data <- LayerData(seurat_subset, layer = layer_name)[targets, , drop = FALSE]
    means <- rowMeans(plot_data)
    means_list <- c(means_list, means)
  }
  
  plot_data_raw <- matrix(means_list, nrow = length(targets), ncol = 2)
  colnames(plot_data_raw) <- cell_types
  rownames(plot_data_raw) <- targets
  
  # 保存到总表
  for (i in 1:nrow(plot_data_raw)) {
    all_csv_data <- rbind(all_csv_data, data.frame(
      miRNA = mirna,
      Target = rownames(plot_data_raw)[i],
      Macrophage = plot_data_raw[i, "Macrophage"],
      Neutrophil = plot_data_raw[i, "Neutrophil"],
      stringsAsFactors = FALSE
    ))
  }
  
  # 过滤低表达基因
  plot_data <- plot_data_raw[rowSums(plot_data_raw) > 0, , drop = FALSE]
  plot_data <- plot_data[apply(plot_data, 1, function(x) any(x >= EXPR_THRESHOLD)), , drop = FALSE]
  if (nrow(plot_data) == 0) next
  
  max_expr_all <- max(max_expr_all, max(plot_data))
  plot_data_list[[mirna]] <- plot_data
}

# 保存总表达量表
csvfilename <- file.path(dirname(RDS_FILE), paste0(sample_name, "_all_target_expression.csv"))
write.csv(all_csv_data, csvfilename, row.names = FALSE)
cat("Saved CSV:", basename(csvfilename), "\n\n")

# ===================== 绘制小提琴图 =====================
cat("=== Generating violin plots ===\n")
for (mirna in target_mirnas) {
  targets <- microrna_df$Symbol[microrna_df$miRNA.ID == mirna]
  targets <- targets[targets %in% rownames(seurat_obj)]
  targets <- targets[tolower(targets) != "actg1"]
  if (length(targets) == 0) {
    cat(mirna, "- no targets found\n")
    next
  }

  plot_data_list_ct <- list()
  for (cell_type in cell_types) {
    seurat_subset <- seurat_obj[, seurat_obj$CellType == cell_type]
    layer_name <- Layers(seurat_subset)[1]
    expr_data <- LayerData(seurat_subset, layer = layer_name)[targets, , drop = FALSE]
    expr_data <- as.data.frame(t(expr_data))
    expr_data$CellType <- cell_type
    plot_data_list_ct[[cell_type]] <- expr_data
  }

  combined_data <- do.call(rbind, plot_data_list_ct)
  combined_data$Cell <- rownames(combined_data)

  plot_df <- data.frame(stringsAsFactors = FALSE)
  for (target in targets) {
    if (!(target %in% colnames(combined_data))) next
    mac_vals <- combined_data[[target]][combined_data$CellType == "Macrophage"]
    neu_vals <- combined_data[[target]][combined_data$CellType == "Neutrophil"]
    all_vals <- c(mac_vals, neu_vals)
    expr_min <- min(all_vals, na.rm = TRUE)
    expr_max <- max(all_vals, na.rm = TRUE)
    if (expr_max > expr_min) {
      mac_norm <- (mac_vals - expr_min) / (expr_max - expr_min)
      neu_norm <- (neu_vals - expr_min) / (expr_max - expr_min)
    } else {
      mac_norm <- mac_vals
      neu_norm <- neu_vals
    }
    mac_df <- data.frame(CellType = "Macrophage", Expression = mac_norm, Target = target, stringsAsFactors = FALSE)
    neu_df <- data.frame(CellType = "Neutrophil", Expression = neu_norm, Target = target, stringsAsFactors = FALSE)
    plot_df <- rbind(plot_df, mac_df, neu_df)
  }

  if (nrow(plot_df) == 0) {
    cat(mirna, "- no valid expression\n")
    next
  }

  plot_df$Target <- factor(plot_df$Target, levels = targets)
  n_targets <- length(unique(plot_df$Target))
  n_cols <- min(3, n_targets)
  n_rows <- ceiling(n_targets / n_cols)

  p <- ggplot(plot_df, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_violin(trim = FALSE) +
    scale_fill_manual(values = c("Macrophage" = "blue", "Neutrophil" = "red")) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(mirna, "Target Expression"), x = "Cell Type", y = "Expression") +
    facet_wrap(~Target, ncol = n_cols, scales = "free_y")

  pngfilename <- file.path(dirname(RDS_FILE), paste0(sample_name, "_", gsub("-", "_", mirna), "_violin.png"))
  png(pngfilename, width = n_cols * 3, height = n_rows * 2.5, units = "in", res = 300)
  print(p)
  dev.off()
  cat("Saved:", basename(pngfilename), "\n")
}

# ===================== 筛选高表达靶基因（修复核心错误） =====================
cat("\n=== Filtering and combining targets ===\n")

# 直接使用已生成的all_csv_data，避免重复计算
all_targets_expr <- all_csv_data
colnames(all_targets_expr)[2] <- "Gene"
all_targets_expr$Target <- paste0(all_targets_expr$miRNA, "|", all_targets_expr$Gene)

cat("Total rows collected:", nrow(all_targets_expr), "\n")

if (nrow(all_targets_expr) == 0) {
  cat("No expression data found. Check miRNA names and gene names.\n")
} else {
  agg_expr <- all_targets_expr
  agg_expr$MaxExpr <- pmax(agg_expr$Macrophage, agg_expr$Neutrophil)
  agg_expr <- agg_expr[agg_expr$MaxExpr >= EXPR_THRESHOLD, ]
  agg_expr <- agg_expr[order(agg_expr$miRNA, agg_expr$MaxExpr, decreasing = TRUE), ]
  cat("Targets with expression >= 0.5:", nrow(agg_expr), "\n\n")

  if (nrow(agg_expr) > 0) {
    cat("=== Generating combined heatmap ===\n")
    heatmap_data <- agg_expr[, c("Target", "Macrophage", "Neutrophil")]
    rownames(heatmap_data) <- heatmap_data$Target
    heatmap_data$Target <- NULL
    heatmap_data <- as.matrix(heatmap_data)

    n_targets <- nrow(heatmap_data)
    height <- max(6, n_targets * 0.4)
    width <- 6

    row_diff <- heatmap_data[, "Macrophage"] - heatmap_data[, "Neutrophil"]
    max_abs <- max(abs(row_diff), na.rm = TRUE)
    break_vals <- seq(-max_abs, max_abs, length.out = 50)

    pngfilename <- file.path(dirname(RDS_FILE), paste0(sample_name, "_combined_targets_heatmap.png"))
    png(pngfilename, width = width, height = height, units = "in", res = 300)
    pheatmap(heatmap_data,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             legend = TRUE,
             breaks = break_vals,
             color = colorRampPalette(c("blue", "white", "red"))(length(break_vals) - 1),
             main = "Target Gene Expression (Macrophage vs Neutrophil)\nMacrophage <- | -> Neutrophil")
    dev.off()
    cat("Saved:", basename(pngfilename), "\n")

    cat("=== Generating bubble plot ===\n")
    bubble_data <- agg_expr[, c("miRNA", "Gene", "Macrophage", "Neutrophil")]
    bubble_data <- bubble_data[!tolower(bubble_data$Gene) %in% c("actg1", "wnt"), ]
    bubble_long <- reshape2::melt(bubble_data, id.vars = c("miRNA", "Gene"),
                                   variable.name = "CellType", value.name = "Expression")
    bubble_norm <- bubble_long %>%
      group_by(Gene) %>%
      mutate(Expression = (Expression - min(Expression, na.rm = TRUE)) /
               (max(Expression, na.rm = TRUE) - min(Expression, na.rm = TRUE)))
    bubble_norm$Label <- paste0(bubble_norm$miRNA, "|", bubble_norm$Gene)
    bubble_norm <- bubble_norm[order(bubble_norm$miRNA, bubble_norm$Gene), ]

    p <- ggplot(bubble_norm, aes(x = Label, y = CellType)) +
      geom_point(aes(color = Expression), alpha = 0.8, size = 4) +
      scale_color_gradient(low = "#f5d375", high = "#FA9FB5") +
      scale_x_discrete(labels = function(x) gsub(".*\\|", "", x)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      labs(title = "Target Expression", x = "Target Gene", y = "Cell Type", color = "Expression")

    svgfilename <- file.path(dirname(RDS_FILE), paste0(sample_name, "_combined_targets_bubble.svg"))
    svglite::svglite(svgfilename, width = max(12, length(unique(bubble_norm$Label)) * 0.5) * 0.8, height = 6)
    print(p)
    dev.off()
    cat("Saved:", basename(svgfilename), "\n")

    bubble_csv <- bubble_norm[, c("miRNA", "Gene", "CellType", "Expression")]
    csvfilename <- file.path(dirname(RDS_FILE), paste0(sample_name, "_combined_targets_bubble_norm.csv"))
    write.csv(bubble_csv, csvfilename, row.names = FALSE)
    cat("Saved:", basename(csvfilename), "\n")

    csvfilename <- file.path(dirname(RDS_FILE), paste0(sample_name, "_combined_targets_heatmap_data.csv"))
    write.csv(agg_expr, csvfilename, row.names = FALSE)
    cat("Saved:", basename(csvfilename), "\n")
  }
}

}

cat("\nDone!\n")