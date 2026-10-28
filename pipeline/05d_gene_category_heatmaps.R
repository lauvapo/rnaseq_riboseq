# ============================================================
# Step 5d: Gene category heatmaps - 5 biological themes
#
# Creates publication-quality paired figures using svglite format
# showing log2FoldChange for RNA-seq and Ribo-seq side-by-side
#
# Outputs:
#   results/figure_themes_1.svg/png (DNA damage + Aromatic AA)
#   results/figure_themes_2.svg/png (Translation + Redox/Lipid/Envelope)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(svglite)
  library(patchwork)
})

# Load all_info
all_info <- read.delim(file.path(dir_sig_genes, "all_info.tsv"))

# Define 5 biological themes
gene_categories <- list(
  list(
    name = "DNAdamage_SOS",
    letter = "A",
    title = "DNA damage / SOS response",
    genes = c(
      "Rv3490", "Rv3395c", "Rv3297", "Rv3296", "Rv3202c", "Rv3201c",
      "Rv3074", "Rv2979c", "Rv2943", "Rv2894c", "Rv2792c", "Rv2791c",
      "Rv2737c", "Rv2594c", "Rv2593c", "Rv2592c", "Rv1000c", "Rv0938",
      "Rv0922", "Rv0605"
    )
  ),
  list(
    name = "Translation_Ribosome",
    letter = "B",
    title = "Translation / ribosome machinery",
    genes = c(
      "Rv0053", "Rv0055", "Rv0056", "Rv0634B", "Rv0651", "Rv0652",
      "Rv0682", "Rv0683", "Rv0684", "Rv0685", "Rv0700", "Rv0701",
      "Rv0702", "Rv0703", "Rv0704", "Rv0705", "Rv0706", "Rv0707",
      "Rv0708", "Rv0709", "Rv0710", "Rv0714", "Rv0715", "Rv0716",
      "Rv0717", "Rv0718", "Rv0719", "Rv0720", "Rv0721", "Rv1014c",
      "Rv1015c", "Rv1641", "Rv1650", "Rv2357c", "Rv2412", "Rv2441c",
      "Rv2442c", "Rv2572c", "Rv2890c", "Rv2904c", "Rv3011c", "Rv3443c",
      "Rv3456c", "Rv3458c", "Rv3459c", "Rv3460c", "Rv3462c", "Rv3907c", "Rv3923c", "Rv3924c"
    )
  ),
  list(
    name = "Redox_Lipid_Envelope",
    letter = "C",
    title = "Redox, lipid & cell-envelope metabolism",
    genes = c(
      "Rv1908c", "Rv2633c", "Rv3846",
      "Rv3250c", "Rv3251c", "Rv3252c",
      "Rv1886c", "Rv3634c", "Rv3794", "Rv3795", "Rv0635", "Rv0636",
      "Rv0637", "Rv0642c", "Rv0643c", "Rv0644c", "Rv3281",
      "Rv0526", "Rv0527", "Rv0528", "Rv2877c"
    )
  ),
  list(
    name = "Aromatic_AA",
    letter = "D",
    title = "Aromatic amino acid biosynthesis",
    genes = c(
      "Rv1603", "Rv1612", "Rv1613", "Rv2538c", "Rv2539c", "Rv2552c",
      "Rv3250c", "Rv3251c", "Rv3252c", "Rv0948c", "Rv1885c"
    )
  )
)

# Function to create individual heatmap plot
create_heatmap_plot <- function(category_info, all_info_df, letter) {

  category_name <- category_info$name
  category_title <- category_info$title
  gene_list <- category_info$genes

  genes_found <- all_info_df %>%
    filter(name %in% gene_list) %>%
    arrange(factor(name, levels = gene_list))

  if (nrow(genes_found) == 0) {
    return(NULL)
  }

  plot_data <- genes_found %>%
    select(gene_symbol, logFC_rna, logFC_ribo) %>%
    pivot_longer(
      cols = starts_with("logFC"),
      names_to = "Assay",
      names_prefix = "logFC_",
      values_to = "logFC"
    ) %>%
    mutate(
      Assay = ifelse(Assay == "rna", "RNA-seq", "Ribo-seq"),
      Assay = factor(Assay, levels = c("RNA-seq", "Ribo-seq")),
      gene_label = factor(gene_symbol, levels = unique(genes_found$gene_symbol))
    ) %>%
    select(-gene_symbol)

  # Create heatmap-style plot
  p <- ggplot(plot_data, aes(x = Assay, y = gene_label, fill = logFC)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#0D88FFFF",
      mid = "white",
      high = "#E16462FF",
      midpoint = 0,
      limits = c(-max(abs(plot_data$logFC), na.rm = TRUE),
                 max(abs(plot_data$logFC), na.rm = TRUE)),
      name = "log2FC"
    ) +
    labs(
      title = paste0(letter, ". ", category_title),
      x = "",
      y = ""
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(
        size = 13, face = "bold", hjust = 0, margin = margin(b = 7)
      ),
      axis.text.x = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 8.5),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.2, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  return(p)
}

# Create all individual plots
cat("[05d_heatmap] Creating paired theme heatmap figures (svglite format)...\n")

plots_list <- list()
for (i in seq_along(gene_categories)) {
  p <- create_heatmap_plot(gene_categories[[i]], all_info, gene_categories[[i]]$letter)
  if (!is.null(p)) {
    plots_list[[gene_categories[[i]]$name]] <- p
  }
}

# Create Figure 1: DNA damage + Aromatic AA (side-by-side)
fig1 <- plots_list$DNAdamage_SOS + plots_list$Aromatic_AA +
  plot_layout(guides = "collect", ncol = 2, widths = c(1, 1)) &
  theme(legend.position = "right")

# Export Figure 1 as PNG
ggsave(
  filename = file.path(dir_results, "figure_themes_1.png"),
  plot = fig1,
  width = 18,
  height = 8.5,
  dpi = 300,
  bg = "white"
)
cat("[05d_heatmap] Figure 1 PNG exported: figure_themes_1.png\n")

# Export Figure 1 as SVG with svglite
svglite::svglite(
  file = file.path(dir_results, "figure_themes_1.svg"),
  width = 18,
  height = 8.5,
  bg = "white"
)
print(fig1)
dev.off()
cat("[05d_heatmap] Figure 1 SVG exported: figure_themes_1.svg\n")

# Create Figure 2: Translation + Redox/Lipid/Envelope (side-by-side)
fig2 <- plots_list$Translation_Ribosome + plots_list$Redox_Lipid_Envelope +
  plot_layout(guides = "collect", ncol = 2, widths = c(1, 1)) &
  theme(legend.position = "right")

# Export Figure 2 as PNG
ggsave(
  filename = file.path(dir_results, "figure_themes_2.png"),
  plot = fig2,
  width = 18,
  height = 16,
  dpi = 300,
  bg = "white"
)
cat("[05d_heatmap] Figure 2 PNG exported: figure_themes_2.png\n")

# Export Figure 2 as SVG with svglite
svglite::svglite(
  file = file.path(dir_results, "figure_themes_2.svg"),
  width = 18,
  height = 16,
  bg = "white"
)
print(fig2)
dev.off()
cat("[05d_heatmap] Figure 2 SVG exported: figure_themes_2.svg\n")

cat("[05d_heatmap] Paired theme heatmap figures complete (svglite format).\n")
