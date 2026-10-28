# ============================================================
# Step 5b: Volcano plots - RNA-seq and Ribo-seq separate + combined
#
# Produces publication-ready volcano plots showing differential
# expression results for RNA-seq and Ribo-seq. Creates:
# - Individual plots for RNA-seq and Ribo-seq
# - Combined faceted plot showing both
# Points colored by significance/direction (upregulated/downregulated/
# not significant). Labeled genes customizable.
#
# Outputs:
#   results/volcano_rna.png/svg
#   results/volcano_ribo.png/svg
#   results/volcano_combined.png/svg
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# Load all_info
all_info <- read.delim(file.path(dir_sig_genes, "all_info.tsv"))

# Prepare RNA-seq data
rna_df <- all_info %>%
  mutate(
    logpadj = -log10(padj_rna + 1e-10),
    Category = case_when(
      padj_rna < fdr_cutoff & logFC_rna > logfc_cutoff ~ "Upregulated",
      padj_rna < fdr_cutoff & logFC_rna < -logfc_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    Dataset = "RNA-seq"
  ) %>%
  select(name, gene_symbol, logFC_rna, padj_rna, logpadj, Category, Dataset)

colnames(rna_df) <- c("gene", "gene_symbol", "log2FoldChange", "padj", "logpadj", "Category", "Dataset")

# Prepare Ribo-seq data
ribo_df <- all_info %>%
  mutate(
    logpadj = -log10(padj_ribo + 1e-10),
    Category = case_when(
      padj_ribo < fdr_cutoff & logFC_ribo > logfc_cutoff ~ "Upregulated",
      padj_ribo < fdr_cutoff & logFC_ribo < -logfc_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    Dataset = "Ribo-seq"
  ) %>%
  select(name, gene_symbol, logFC_ribo, padj_ribo, logpadj, Category, Dataset)

colnames(ribo_df) <- c("gene", "gene_symbol", "log2FoldChange", "padj", "logpadj", "Category", "Dataset")

# Combine datasets
combined_df <- bind_rows(rna_df, ribo_df)
combined_df$Dataset <- factor(combined_df$Dataset, levels = c("RNA-seq", "Ribo-seq"))

# Define genes to label
genes_to_label_all <- c(
  "leuD", "leuC", "alkB",
  "rubA", "rubB",
  "rpsN1", "rpsH", "rpsS", "rpsQ",
  "efpA"
)

genes_to_label <- list(
  "RNA-seq" = genes_to_label_all,
  "Ribo-seq" = genes_to_label_all
)

# ============================================================
# COMBINED VOLCANO PLOT (faceted)
# ============================================================

# Calculate x-axis limits based on plotted genes (those with p.adjust values) with 10% padding
plotted_genes <- combined_df[!is.na(combined_df$padj), ]$log2FoldChange
x_range <- range(plotted_genes, na.rm = TRUE)
x_padding <- diff(x_range) * 0.1
x_limits <- c(x_range[1] - x_padding, x_range[2] + x_padding)

# Prepare label dataframe for combined plot
label_df_combined <- combined_df %>%
  filter(
    (Dataset == "RNA-seq" & gene_symbol %in% genes_to_label[["RNA-seq"]]) |
    (Dataset == "Ribo-seq" & gene_symbol %in% genes_to_label[["Ribo-seq"]])
  )

# Create combined volcano plot
volcano_combined <- ggplot(combined_df, aes(x = log2FoldChange, y = logpadj, color = Category)) +

  geom_point(alpha = 0.7, size = 5) +

  geom_hline(yintercept = -log10(fdr_cutoff), color = "darkgray", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), color = "darkgray", linetype = "dashed", linewidth = 0.8) +

  geom_text_repel(
    data = label_df_combined,
    aes(label = gene_symbol),
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    nudge_y = 0.7,
    max.overlaps = 20
  ) +

  facet_wrap(~Dataset, nrow = 1) +

  scale_color_manual(
    values = c(
      "Not significant" = "gray",
      "Upregulated" = "#FF6961FF",
      "Downregulated" = "#4285F4FF"
    )
  ) +

  scale_x_continuous(limits = x_limits, breaks = pretty(x_limits, 4)) +
  scale_y_continuous(breaks = pretty(combined_df$logpadj, 6)) +

  theme_minimal(base_size = 24) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 30, face = "bold"),
    axis.title = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 25),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30, face = "bold"),
    legend.key.size = unit(1.2, "cm")
  ) +

  labs(
    x = "log2FoldChange",
    y = "-log10(padj)",
    color = "Category"
  )

# Export combined volcano
ggsave(
  filename = file.path(dir_results, "volcano_combined.png"),
  plot = volcano_combined,
  width = 25, height = 15, dpi = 300
)
cat("[05b_volcano] Combined PNG exported: volcano_combined.png\n")

ggsave(
  filename = file.path(dir_results, "volcano_combined.svg"),
  plot = volcano_combined,
  device = "svg",
  width = 10, height = 8, units = "in"
)
cat("[05b_volcano] Combined SVG exported: volcano_combined.svg\n")

# ============================================================
# INDIVIDUAL RNA-SEQ VOLCANO PLOT
# ============================================================

rna_only_df <- filter(rna_df, Dataset == "RNA-seq")

# Calculate x-axis limits based on plotted genes (those with p.adjust values) with 10% padding
rna_plotted <- rna_only_df[!is.na(rna_only_df$padj), ]$log2FoldChange
rna_x_range <- range(rna_plotted, na.rm = TRUE)
rna_x_padding <- diff(rna_x_range) * 0.1
rna_x_limits <- c(rna_x_range[1] - rna_x_padding, rna_x_range[2] + rna_x_padding)

label_df_rna <- filter(rna_only_df, gene_symbol %in% genes_to_label_all)

volcano_rna <- ggplot(rna_only_df, aes(x = log2FoldChange, y = logpadj, color = Category)) +

  geom_point(alpha = 0.7, size = 5) +

  geom_hline(yintercept = -log10(fdr_cutoff), color = "darkgray", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), color = "darkgray", linetype = "dashed", linewidth = 0.8) +

  geom_text_repel(
    data = label_df_rna,
    aes(label = gene_symbol),
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    nudge_y = 0.7,
    max.overlaps = 20
  ) +

  scale_color_manual(
    values = c(
      "Not significant" = "gray",
      "Upregulated" = "#FF6961FF",
      "Downregulated" = "#4285F4FF"
    )
  ) +

  scale_x_continuous(limits = rna_x_limits, breaks = pretty(rna_x_limits, 6)) +
  scale_y_continuous(breaks = pretty(rna_only_df$logpadj, 6)) +

  theme_minimal(base_size = 24) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 25),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30, face = "bold"),
    legend.key.size = unit(1.2, "cm")
  ) +

  labs(
    title = "RNA-seq Volcano Plot",
    x = "log2FoldChange",
    y = "-log10(padj)",
    color = "Category"
  )

# Export RNA-seq volcano
ggsave(
  filename = file.path(dir_results, "volcano_rna.png"),
  plot = volcano_rna,
  width = 14, height = 12, dpi = 300
)
cat("[05b_volcano] RNA-seq PNG exported: volcano_rna.png\n")

ggsave(
  filename = file.path(dir_results, "volcano_rna.svg"),
  plot = volcano_rna,
  device = "svg",
  width = 11, height = 9, units = "in"
)
cat("[05b_volcano] RNA-seq SVG exported: volcano_rna.svg\n")

# ============================================================
# INDIVIDUAL RIBO-SEQ VOLCANO PLOT
# ============================================================

ribo_only_df <- filter(ribo_df, Dataset == "Ribo-seq")

# Calculate x-axis limits based on plotted genes (those with p.adjust values) with 10% padding
ribo_plotted <- ribo_only_df[!is.na(ribo_only_df$padj), ]$log2FoldChange
ribo_x_range <- range(ribo_plotted, na.rm = TRUE)
ribo_x_padding <- diff(ribo_x_range) * 0.1
ribo_x_limits <- c(ribo_x_range[1] - ribo_x_padding, ribo_x_range[2] + ribo_x_padding)

label_df_ribo <- filter(ribo_only_df, gene_symbol %in% genes_to_label_all)

volcano_ribo <- ggplot(ribo_only_df, aes(x = log2FoldChange, y = logpadj, color = Category)) +

  geom_point(alpha = 0.7, size = 5) +

  geom_hline(yintercept = -log10(fdr_cutoff), color = "darkgray", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), color = "darkgray", linetype = "dashed", linewidth = 0.8) +

  geom_text_repel(
    data = label_df_ribo,
    aes(label = gene_symbol),
    size = 10,
    fontface = "bold",
    show.legend = FALSE,
    nudge_y = 0.7,
    max.overlaps = 20
  ) +

  scale_color_manual(
    values = c(
      "Not significant" = "gray",
      "Upregulated" = "#FF6961FF",
      "Downregulated" = "#4285F4FF"
    )
  ) +

  scale_x_continuous(limits = ribo_x_limits, breaks = pretty(ribo_x_limits, 6)) +
  scale_y_continuous(breaks = pretty(ribo_only_df$logpadj, 6)) +

  theme_minimal(base_size = 24) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 25),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30, face = "bold"),
    legend.key.size = unit(1.2, "cm")
  ) +

  labs(
    title = "Ribo-seq Volcano Plot",
    x = "log2FoldChange",
    y = "-log10(padj)",
    color = "Category"
  )

# Export Ribo-seq volcano
ggsave(
  filename = file.path(dir_results, "volcano_ribo.png"),
  plot = volcano_ribo,
  width = 14, height = 12, dpi = 300
)
cat("[05b_volcano] Ribo-seq PNG exported: volcano_ribo.png\n")

ggsave(
  filename = file.path(dir_results, "volcano_ribo.svg"),
  plot = volcano_ribo,
  device = "svg",
  width = 11, height = 9, units = "in"
)
cat("[05b_volcano] Ribo-seq SVG exported: volcano_ribo.svg\n")
