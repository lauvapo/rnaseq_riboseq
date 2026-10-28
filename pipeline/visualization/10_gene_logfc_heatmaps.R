# ============================================================
# Step 10: Gene-level RNA-seq vs Ribo-seq heatmap
#
# Steps 13/14 heatmap GO/pathway TERMS. This step goes one level below
# that, to individual GENES: for each cross-assay DE category
# (Co-regulated, Only RNA-seq, Only Ribo-seq), the top genes by
# |log2FC|, sorted by fold-change (not clustered, not alphabetical),
# RNA-seq and Ribo-seq as the two columns, diverging blue-red gradient.
# The direct answer to "how do these two assays compare for the genes
# that moved the most" - a complement to summary_transcriptome_vs_
# translatome.png (same comparison, scatter form) and
# summary_holistic_gene_view.png (same comparison, hand-picked genes,
# bar form).
#
# 696 genes are DE-significant in total (see summary_de_categories.png)
# - too many for one legible heatmap - so this caps each category's
# panel at TOP_N_GENES by effect size rather than showing everything.
#
# Reuses theme_pub()/save_pub() from 12_summary_figures.R, which always
# runs before this step in run_pipeline.R.
#
# Outputs (results/):
#   gene_logfc_heatmap.{png,pdf}   3 panels (A: Co-regulated, B: Only
#                                   RNA-seq, C: Only Ribo-seq)
# ============================================================

all_info_gene_hm <- read.delim(file.path(dir_sig_genes, "all_info.tsv"))

# Shows ALL genes in this category (using whichever assay has a value -
# RNA-seq preferred, falling back to Ribo-seq for genes where only
# Ribo-seq is significant/available), ordered by that same signed
# fold-change so the panel reads as one continuous ramp from most up-
# to most down-regulated, with color reinforcing direction.
make_gene_panel <- function(category_name, limit) {
  d <- all_info_gene_hm %>%
    filter(category == category_name, !is.na(gene_symbol)) %>%
    mutate(sort_fc = ifelse(!is.na(logFC_rna), logFC_rna, logFC_ribo)) %>%
    filter(!is.na(sort_fc)) %>%
    arrange(desc(sort_fc))

  if (nrow(d) == 0) return(NULL)

  # name (locus tag) is the unique key throughout this project;
  # gene_symbol is display-only and can collide for genes with no
  # assigned symbol (both fall back to their locus tag) - disambiguate
  # any collision rather than silently merging two different genes into
  # one heatmap row.
  d$row_label <- ifelse(duplicated(d$gene_symbol) | duplicated(d$gene_symbol, fromLast = TRUE),
                         paste0(d$gene_symbol, " (", d$name, ")"), d$gene_symbol)
  row_order <- rev(d$row_label)

  plot_df <- d %>%
    select(row_label, logFC_rna, logFC_ribo) %>%
    tidyr::pivot_longer(c(logFC_rna, logFC_ribo), names_to = "assay", values_to = "logFC") %>%
    mutate(
      assay = factor(ifelse(assay == "logFC_rna", "RNA-seq", "Ribo-seq"), levels = c("RNA-seq", "Ribo-seq")),
      row_label = factor(row_label, levels = row_order)
    )

  ggplot(plot_df, aes(x = assay, y = row_label, fill = logFC)) +
    geom_tile(color = "white", linewidth = 0.5, na.rm = FALSE) +
    scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#993C1D", midpoint = 0,
                          limits = c(-limit, limit), oob = scales::squish, na.value = "grey88",
                          name = "log2FC") +
    labs(subtitle = paste0(category_name, " (top ", nrow(d), " by |log2FC|)"), x = NULL, y = NULL) +
    theme_pub(21) +
    theme(panel.grid = element_blank(),
          plot.subtitle = element_text(face = "bold", color = "black", size = rel(1.05)),
          axis.text.y = element_text(face = "italic", size = rel(0.85)))
}

category_order <- c("Co-regulated", "Only RNA-Seq", "Only Ribo-Seq")

limit_gene_hm <- all_info_gene_hm %>%
  filter(category %in% category_order) %>%
  summarise(m = max(abs(c(logFC_rna, logFC_ribo)), na.rm = TRUE)) %>%
  pull(m)

panels_gene_hm <- lapply(category_order, make_gene_panel, limit = limit_gene_hm)
keep_gene_hm <- !vapply(panels_gene_hm, is.null, logical(1))

if (any(keep_gene_hm)) {
  panels_gene_hm <- panels_gene_hm[keep_gene_hm]
  heights_gene_hm <- vapply(panels_gene_hm, function(p) length(levels(p$data$row_label)), integer(1)) + 3.0

  p_gene_hm <- wrap_plots(panels_gene_hm, ncol = 1, heights = heights_gene_hm) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "RNA-seq vs Ribo-seq: top genes by fold-change",
      subtitle = "Genes ranked by |log2FC|, ordered by signed fold-change",
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 25),
                    plot.subtitle = element_text(color = "grey40", size = 17),
                    plot.tag = element_text(face = "bold", size = 21))
    )

  # Dynamically size the figure based on number of genes
  # Use 0.08 inches per gene row for all genes
  dynamic_height <- min(50, max(8, 3.0 + 0.08 * sum(heights_gene_hm)))
  save_pub(p_gene_hm, "gene_logfc_heatmap", width = 8.5, height = dynamic_height)
} else {
  message("[15_gene_logfc_heatmaps] No genes found in any category - skipping.")
}

cat("[15_gene_logfc_heatmaps] Done.\n")
