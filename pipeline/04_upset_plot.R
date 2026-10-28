# ============================================================
# Step 4: UpSet plot - gene set intersections
#
# Creates an UpSet plot showing the intersections between:
#   - RNA-seq Up-regulated genes
#   - RNA-seq Down-regulated genes
#   - Ribo-seq Up-regulated genes
#   - Ribo-seq Down-regulated genes
#
# Uses ComplexUpset with customized styling for publication quality
#
# Outputs:
#   results/upset_rna_ribo.png (1400 x 700 px)
#   results/upset_rna_ribo.pdf (14 x 7 in)
# ============================================================

suppressPackageStartupMessages({
  library(ComplexUpset)
  library(ggplot2)
  library(dplyr)
})

# Load all_info which contains DE results for both assays
all_info <- read.delim(file.path(dir_sig_genes, "all_info.tsv"))

# Define gene sets by direction and significance
rna_up <- all_info %>%
  filter(!is.na(padj_rna), padj_rna < fdr_cutoff, logFC_rna > logfc_cutoff) %>%
  pull(name)

rna_down <- all_info %>%
  filter(!is.na(padj_rna), padj_rna < fdr_cutoff, logFC_rna < -logfc_cutoff) %>%
  pull(name)

ribo_up <- all_info %>%
  filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo > logfc_cutoff) %>%
  pull(name)

ribo_down <- all_info %>%
  filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo < -logfc_cutoff) %>%
  pull(name)

# Create a binary matrix for ComplexUpset
genes <- unique(c(rna_up, rna_down, ribo_up, ribo_down))
upset_df <- data.frame(
  gene = genes,
  check.names = FALSE,
  "RNA-seq Up" = genes %in% rna_up,
  "RNA-seq Down" = genes %in% rna_down,
  "Ribo-seq Up" = genes %in% ribo_up,
  "Ribo-seq Down" = genes %in% ribo_down,
  row.names = genes
)

# Create gene list for reference
gene_list <- list(
  "RNA-seq Up" = rna_up,
  "RNA-seq Down" = rna_down,
  "Ribo-seq Up" = ribo_up,
  "Ribo-seq Down" = ribo_down
)

# Create the upset plot with custom styling
upset_plot <- ComplexUpset::upset(
  upset_df,
  intersect = names(gene_list),
  min_size = 5,
  matrix = ComplexUpset::intersection_matrix(
    geom = geom_point(size = 8),
    segment = geom_segment(linewidth = 4)
  ),
  set_sizes = ComplexUpset::upset_set_size() +
    theme(
      axis.text.x = element_text(size = 22),
      text = element_text(size = 30, face = "bold")
    ),
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(
      width = 0.6,
      text = list(
        size = 13,
        fontface = "bold"
      ),
      mapping = aes(fill = 'bars_color')
    ) +
    scale_fill_manual(
      values = c('bars_color' = 'black'),
      guide = 'none'
    ) +
    theme(
      axis.text.y = element_text(size = 24, face = 'plain'),
      text = element_text(size = 30, face = 'bold')
    )
  ),
  sort_sets = FALSE,
  sort_intersections_by = "cardinality",
  queries = list(
    ComplexUpset::upset_query(set = "RNA-seq Up", fill = "#E16462FF", color = "#E16462FF"),
    ComplexUpset::upset_query(set = "RNA-seq Down", fill = "#0D0887FF", color = "#0D0887FF"),
    ComplexUpset::upset_query(set = "Ribo-seq Up", fill = "#E16462FF", color = "#E16462FF"),
    ComplexUpset::upset_query(set = "Ribo-seq Down", fill = "#0D0887FF", color = "#0D0887FF")
  )
) +
theme(
  text = element_blank(),
  axis.text.y = element_text(size = 35, face = "bold", color = 'black')
)

# Export PNG
png(file.path(dir_results, "upset_rna_ribo.png"),
    width = 1400, height = 700)
print(upset_plot)
dev.off()
cat("[04_upset] PNG exported: upset_rna_ribo.png\n")

# Export SVG
svg(file.path(dir_results, "upset_rna_ribo.svg"),
    width = 14, height = 7)
print(upset_plot)
dev.off()
cat("[04_upset] SVG exported: upset_rna_ribo.svg\n")

# Export PDF
pdf(file.path(dir_results, "upset_rna_ribo.pdf"),
    width = 14, height = 7)
print(upset_plot)
dev.off()
cat("[04_upset] PDF exported: upset_rna_ribo.pdf\n")

# Print summary
cat("[04_upset] Gene set sizes:\n")
cat("  RNA-seq Up:    ", length(rna_up), "\n")
cat("  RNA-seq Down:  ", length(rna_down), "\n")
cat("  Ribo-seq Up:   ", length(ribo_up), "\n")
cat("  Ribo-seq Down: ", length(ribo_down), "\n")
