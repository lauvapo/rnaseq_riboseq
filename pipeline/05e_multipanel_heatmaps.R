# ============================================================
# Step 5e: Multipanel heatmap figure
#
# Combines all 11 curated gene category heatmaps into a single
# publication-quality multipanel figure with legible gene names
#
# Outputs:
#   results/figure_all_categories_heatmaps.png
#   results/figure_all_categories_heatmaps.svg
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

# Define all gene categories
gene_categories <- list(
  list(
    name = "SOS_DNA",
    title = "A. SOS Response & DNA Repair",
    genes = c(
      "Rv3490", "Rv3395c", "Rv3297", "Rv3296", "Rv3202c", "Rv3201c",
      "Rv3074", "Rv2979c", "Rv2943", "Rv2894c", "Rv2792c", "Rv2791c",
      "Rv2737c", "Rv2594c", "Rv2593c", "Rv2592c", "Rv1000c", "Rv0938",
      "Rv0922", "Rv0605"
    )
  ),
  list(
    name = "Ribosome",
    title = "B. Ribosomal structure",
    genes = c(
      "Rv0053", "Rv0055", "Rv0056", "Rv0634B", "Rv0651", "Rv0652",
      "Rv0682", "Rv0683", "Rv0700", "Rv0701", "Rv0702", "Rv0703",
      "Rv0704", "Rv0705", "Rv0706", "Rv0707", "Rv0708", "Rv0709",
      "Rv0710", "Rv0714", "Rv0715", "Rv0716", "Rv0717", "Rv0718",
      "Rv0719", "Rv0720", "Rv0721", "Rv1015c", "Rv2412", "Rv2441c",
      "Rv2442c", "Rv2890c", "Rv2904c", "Rv3443c", "Rv3456c", "Rv3458c",
      "Rv3459c", "Rv3460c", "Rv3924c"
    )
  ),
  list(
    name = "Translation",
    title = "C. Translation machinery",
    genes = c(
      "Rv0053", "Rv0055", "Rv0056", "Rv0634B", "Rv0651", "Rv0652",
      "Rv0682", "Rv0683", "Rv0684", "Rv0685", "Rv0700", "Rv0701",
      "Rv0702", "Rv0703", "Rv0704", "Rv0705", "Rv0706", "Rv0707",
      "Rv0708", "Rv0709", "Rv0710", "Rv0714", "Rv0715", "Rv0716",
      "Rv0717", "Rv0718", "Rv0719", "Rv0720", "Rv0721", "Rv1014c",
      "Rv1015c", "Rv1641", "Rv1650", "Rv2357c", "Rv2412", "Rv2441c",
      "Rv2442c", "Rv2572c", "Rv2890c", "Rv2904c", "Rv3011c", "Rv3443c",
      "Rv3456c", "Rv3458c", "Rv3459c", "Rv3460c", "Rv3462c", "Rv3924c"
    )
  ),
  list(
    name = "tRNA",
    title = "D. tRNA binding",
    genes = c(
      "Rv0682", "Rv0683", "Rv0700", "Rv0708", "Rv0716", "Rv1014c",
      "Rv1650", "Rv3460c", "Rv3907c", "Rv3923c"
    )
  ),
  list(
    name = "Aromatic_AA",
    title = "E. Aromatic amino acid biosynthesis",
    genes = c(
      "Rv1603", "Rv1612", "Rv1613", "Rv2538c", "Rv2539c", "Rv2552c",
      "Rv3250c", "Rv3251c", "Rv3252c", "Rv0948c", "Rv1885c"
    )
  ),
  list(
    name = "Mycolate",
    title = "F. Mycolate biosynthesis",
    genes = c(
      "Rv1886c", "Rv3634c", "Rv3794", "Rv3795", "Rv0635", "Rv0636",
      "Rv0637", "Rv0642c", "Rv0643c", "Rv0644c"
    )
  ),
  list(
    name = "CellEnvelope",
    title = "G. Cell envelope biosynthesis",
    genes = c(
      "Rv1886c", "Rv3634c", "Rv3794", "Rv3795", "Rv0635", "Rv0636",
      "Rv0637", "Rv0642c", "Rv0643c", "Rv0644c", "Rv3281"
    )
  ),
  list(
    name = "Redox",
    title = "H. Redox & antioxidant defense",
    genes = c(
      "Rv1908c", "Rv2633c", "Rv3846"
    )
  ),
  list(
    name = "Lipid",
    title = "I. Fatty acid & lipid metabolism",
    genes = c(
      "Rv3250c", "Rv3251c", "Rv3252c"
    )
  ),
  list(
    name = "ElectronTransport",
    title = "J. Electron transport & respiration",
    genes = c(
      "Rv0526", "Rv0527", "Rv0528", "Rv2877c"
    )
  )
)

# Color palette
col_fun_gradient <- function(logfc_vals) {
  colorRamp2(
    c(-max(abs(logfc_vals), na.rm = TRUE),
      0,
      max(abs(logfc_vals), na.rm = TRUE)),
    c("#0D88FFFF", "white", "#E16462FF")
  )
}

# Function to create individual heatmap plot
create_heatmap_plot <- function(category_info, all_info_df, include_legend = FALSE) {

  category_name <- category_info$name
  category_title <- category_info$title
  gene_list <- category_info$genes

  genes_found <- all_info_df %>%
    filter(name %in% gene_list) %>%
    arrange(factor(name, levels = gene_list))

  if (nrow(genes_found) == 0) return(NULL)

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

  p <- ggplot(plot_data, aes(x = Assay, y = gene_label, fill = logFC)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "#0D88FFFF",
      mid = "white",
      high = "#E16462FF",
      midpoint = 0,
      limits = c(-max(abs(plot_data$logFC), na.rm = TRUE),
                 max(abs(plot_data$logFC), na.rm = TRUE)),
      name = "log2FC",
      guide = if (include_legend) "colorbar" else "none"
    ) +
    labs(title = category_title, x = "", y = "") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(
        size = 11, face = "bold", hjust = 0, margin = margin(b = 5)
      ),
      axis.text.x = element_text(size = 9, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = if (include_legend) "right" else "none",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.3, "cm")
    )

  return(p)
}

# Create all heatmaps
cat("[05e_multipanel] Creating multipanel heatmap figure...\n")

plots_list <- list()
for (i in seq_along(gene_categories)) {
  cat("[05e_multipanel] Creating panel", LETTERS[i], "...\n")
  include_legend <- (i == length(gene_categories))  # legend only on last panel
  p <- create_heatmap_plot(gene_categories[[i]], all_info, include_legend = include_legend)
  if (!is.null(p)) {
    plots_list[[i]] <- p
  }
}

# Combine into multipanel (4 columns × 3 rows)
multipanel <- wrap_plots(plots_list, ncol = 4, nrow = 3) +
  plot_layout(
    guides = "collect",
    axis_titles = "collect"
  ) &
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

# Export PNG
ggsave(
  filename = file.path(dir_results, "figure_all_categories_heatmaps.png"),
  plot = multipanel,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white"
)
cat("[05e_multipanel] PNG exported: figure_all_categories_heatmaps.png\n")

# Export SVG with svglite
svglite::svglite(
  file = file.path(dir_results, "figure_all_categories_heatmaps.svg"),
  width = 16,
  height = 12,
  bg = "white"
)
print(multipanel)
dev.off()
cat("[05e_multipanel] SVG exported: figure_all_categories_heatmaps.svg\n")

cat("[05e_multipanel] Multipanel heatmap figure complete.\n")
