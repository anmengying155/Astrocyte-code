# ============================================================================
# miRNA Target Score Visualization
# Purpose: Plot miRNA-target scores and export top genes per miRNA
# Input: Top5_RNA_Targets.xlsx
# Output: target_score_scatter.svg, miRNA_target_predictions.csv
# ============================================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(svglite)
library(ggrepel)

# Define paths relative to project root
SCRIPT_DIR <- dirname(parent.frame(2)$ofile)
PROJECT_ROOT <- ifelse(is.null(SCRIPT_DIR), ".", dirname(SCRIPT_DIR))

TARGET_FILE <- file.path(PROJECT_ROOT, "..", "Data", "Top5_RNA_Targets.xlsx")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "..", "Results")

data <- read_excel(TARGET_FILE)
colnames(data)[6] <- "TargetScore"
data <- data[!is.na(data$`TargetScore`), ]
data <- data[!is.na(data$`miRNA ID`), ]
data <- data[data$`miRNA ID` != "NA", ]
cat("Columns:", colnames(data), "\n")
cat("Rows:", nrow(data), "\n")

top_genes <- data %>%
  group_by(`miRNA ID`) %>%
  slice_max(`TargetScore`, n = 3) %>%
  ungroup()

p <- ggplot(data, aes(x = `TargetScore`, y = `miRNA ID`)) +
  geom_point(aes(color = ifelse(`TargetScore` >= 80, "red", "gray50")),
             size = 3) +
  scale_color_manual(values = c("red" = "red", "gray50" = "gray50"), guide = "none") +
  ggrepel::geom_text_repel(
    data = top_genes,
    aes(label = Symbol),
    size = 3,
    max.overlaps = Inf,
    segment.color = NA,
    box.padding = 0.3,
    nudge_x = 1,
    direction = "y"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10)
  ) +
  labs(title = "miRNA-Target Score Plot", x = "TargetScore", y = "miRNA ID")

svgfilename <- file.path(OUTPUT_DIR, "target_score_scatter.svg")
svglite::svglite(svgfilename, width = 12, height = max(8, length(unique(data$`miRNA ID`)) * 0.8))
print(p)
dev.off()
cat("Saved:", basename(svgfilename), "\n")

target_csvfile <- file.path(PROJECT_ROOT, "..", "Data", "miRNA_target_predictions.csv")
write.csv(top_genes, target_csvfile, row.names = FALSE)
cat("Saved:", basename(target_csvfile), "\n")