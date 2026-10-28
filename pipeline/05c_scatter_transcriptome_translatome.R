# ============================================================
# Step 5c: Scatter plot - Transcriptome vs Translatome
#
# Creates a scatter plot comparing RNA-seq and Ribo-seq fold-changes,
# showing co-regulated, opposite, and assay-specific genes with color-coding.
# Optional gene labels for specific genes of interest.
#
# Outputs:
#   results/scatter_translatome.png
#   results/scatter_translatome.svg
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(scales)
})

# Load all_info
all_info <- read.delim(file.path(dir_sig_genes, "all_info.tsv"))

# Prepare data
merged_df <- all_info %>%
  mutate(
    logFC_rna = logFC_rna,
    logFC_ribo = logFC_ribo,
    padj_rna = padj_rna,
    padj_ribo = padj_ribo
  ) %>%
  select(name, gene_symbol, logFC_rna, logFC_ribo, padj_rna, padj_ribo, category)

colnames(merged_df)[1:2] <- c("gene", "gene_symbol")

# Drop genes with padj NA in BOTH assays - DESeq2's independent
# filtering excluded them from significance testing entirely in both
# RNA-seq and Ribo-seq (almost always near-zero raw counts), so their
# logFC is sampling noise rather than signal. Left in, a handful of
# these wild outliers blow out the axis range and compress every real
# point into a tiny corner of the plot.
merged_df <- merged_df %>% filter(!(is.na(padj_rna) & is.na(padj_ribo)))

# Verify category labels exist, if not recreate them
if (!("category" %in% colnames(merged_df)) || is.na(merged_df$category[1])) {
  sig_rna <- !is.na(merged_df$padj_rna) &
    abs(merged_df$logFC_rna) > logfc_cutoff &
    merged_df$padj_rna < fdr_cutoff

  sig_ribo <- !is.na(merged_df$padj_ribo) &
    abs(merged_df$logFC_ribo) > logfc_cutoff &
    merged_df$padj_ribo < fdr_cutoff

  same_direction <- sign(merged_df$logFC_rna) == sign(merged_df$logFC_ribo)
  opposite <- sign(merged_df$logFC_rna) != sign(merged_df$logFC_ribo)

  merged_df$category <- "Not significant"
  merged_df$category[sig_rna & sig_ribo & same_direction] <- "Co-regulated"
  merged_df$category[sig_rna & sig_ribo & opposite] <- "Opposite"
  merged_df$category[!sig_rna & sig_ribo] <- "Only Ribo-Seq"
  merged_df$category[sig_rna & !sig_ribo] <- "Only RNA-Seq"
}

# Define colors using viridis plasma
cols <- viridis_pal(option = "plasma")(9)
color_map <- c(
  "Co-regulated" = cols[5],
  "Only Ribo-Seq" = cols[6],
  "Only RNA-Seq" = cols[3],
  "Opposite" = cols[4],
  "Not significant" = "gray"
)

# Define genes to label
genes_to_label <- c(
  "leuD", "leuC", "alkB",
  "rubA", "rubB",
  "rpsN1", "rpsH", "rpsS", "rpsQ",
  "efpA"
)

# Prepare label dataframe
label_df <- merged_df[merged_df$gene_symbol %in% genes_to_label, ]

# Create scatter plot
scatter_translatome <- ggplot(merged_df, aes(x = logFC_rna, y = logFC_ribo, color = category)) +

  geom_point(alpha = 0.7, size = 7) +

  geom_text_repel(
    data = label_df,
    aes(label = gene_symbol),
    size = 12,
    max.overlaps = 20,
    fontface = "bold",
    show.legend = FALSE
  ) +

  geom_hline(yintercept = c(-logfc_cutoff, logfc_cutoff), color = "gray50", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), color = "gray50", alpha = 0.5, linetype = "dashed") +

  scale_x_continuous(
    breaks = seq(
      floor(min(merged_df$logFC_rna, na.rm = TRUE)),
      ceiling(max(merged_df$logFC_rna, na.rm = TRUE)),
      by = 0.5
    )
  ) +
  scale_y_continuous(
    breaks = seq(
      floor(min(merged_df$logFC_ribo, na.rm = TRUE)),
      ceiling(max(merged_df$logFC_ribo, na.rm = TRUE)),
      by = 0.5
    )
  ) +

  scale_color_manual(
    values = color_map,
    drop = TRUE,
    na.translate = FALSE,
    breaks = c("Co-regulated", "Only Ribo-Seq", "Only RNA-Seq", "Not significant")
  ) +

  theme_minimal(base_size = 24) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 35, face = "bold"),
    axis.title.y = element_text(size = 35, face = "bold"),
    axis.text.x = element_text(size = 35),
    axis.text.y = element_text(size = 35),
    legend.text = element_text(size = 35),
    legend.title = element_text(size = 35, face = "bold"),
    legend.key.size = unit(1.2, "cm"),
    legend.position = "right"
  ) +

  labs(
    x = "log2FoldChange (RNA-seq)",
    y = "log2FoldChange (Ribo-seq)",
    color = "Category"
  )

# Export PNG (high resolution for presentations)
ggsave(
  filename = file.path(dir_results, "scatter_translatome.png"),
  plot = scatter_translatome,
  width = 25, height = 15, dpi = 300
)
cat("[05c_scatter] PNG exported: scatter_translatome.png\n")

# Export SVG (vector for manuscripts)
ggsave(
  filename = file.path(dir_results, "scatter_translatome.svg"),
  plot = scatter_translatome,
  device = "svg",
  width = 20, height = 10, units = "in"
)
cat("[05c_scatter] SVG exported: scatter_translatome.svg\n")

# Print summary
cat("[05c_scatter] Category summary:\n")
print(table(merged_df$category))
