# ============================================================
# Step 9: Term-level log2FC heatmaps
#
# Step 13's heatmaps color each enriched term by how SIGNIFICANT it is
# (-log10 adjusted p). This step asks a different question, matching
# the style of this project's earlier, pre-pipeline exploratory
# heatmaps (../heatmaps/selected/GO_*.png, PW_*.png): for the terms
# significant in a given gene set, what's the actual EFFECT (mean
# log2FC of the term's member genes, in both assays) - i.e. direction
# and magnitude, not confidence.
#
# One figure per gene set: rows are that gene set's significant terms
# (all 4 annotation categories - Biological Process, Cellular
# Component, Molecular Function, Pathway - as separate panels), columns
# are RNA-seq and Ribo-seq, colored by each term's mean log2FC across
# its member genes (diverging scale). Unlike the earlier exploratory
# version, every category is broken out explicitly rather than only
# whichever happened to have hits.
#
# Reuses theme_pub() (from 12) and save_enrichment_fig()/
# annotation_order (from 13), which always run immediately before this
# step in run_pipeline.R.
#
# Outputs (enrichment_results/):
#   logfc_heatmap_go_pathway_{RNAup,RNAdown,RIBOup,RIBOdown}.{png,pdf}
# ============================================================

all_info_logfc <- read.delim(file.path(dir_sig_genes, "all_info.tsv"))

# Adds each term's mean log2FC (across its member genes, from the
# "/"-separated geneID column enricher() writes) in both assays - the
# per-term effect-direction summary this heatmap actually plots.
add_term_logfc <- function(df) {
  gene_lists <- strsplit(df$geneID, "/")
  df$logFC_rna  <- vapply(gene_lists, function(g)
    mean(all_info_logfc$logFC_rna[match(g, all_info_logfc$name)], na.rm = TRUE), numeric(1))
  df$logFC_ribo <- vapply(gene_lists, function(g)
    mean(all_info_logfc$logFC_ribo[match(g, all_info_logfc$name)], na.rm = TRUE), numeric(1))
  df
}

# Builds one panel (one annotation category) for a single gene set:
# terms kept in their original (most-significant-first) order from the
# ORA table, rather than re-sorted or re-clustered - this is a
# complement to step 13's clustered significance view, not a
# replacement, so it deliberately reads differently.
make_logfc_panel <- function(df, annotation_name, limit) {
  d <- df %>% filter(Annotation == annotation_name) %>% distinct(Description, .keep_all = TRUE)
  if (nrow(d) == 0) return(NULL)

  term_levels <- rev(d$Description)
  plot_df <- d %>%
    select(Description, logFC_rna, logFC_ribo) %>%
    tidyr::pivot_longer(-Description, names_to = "assay", values_to = "logFC") %>%
    mutate(
      assay = factor(ifelse(assay == "logFC_rna", "RNA-seq", "Ribo-seq"), levels = c("RNA-seq", "Ribo-seq")),
      label_wrapped = vapply(as.character(Description), function(s) paste(strwrap(s, width = 34), collapse = "\n"), character(1))
    )
  wrapped_levels <- vapply(term_levels, function(s) paste(strwrap(s, width = 34), collapse = "\n"), character(1))
  plot_df$label_wrapped <- factor(plot_df$label_wrapped, levels = wrapped_levels)

  ggplot(plot_df, aes(x = assay, y = label_wrapped, fill = logFC)) +
    geom_tile(color = "white", linewidth = 0.6) +
    scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#993C1D", midpoint = 0,
                          limits = c(-limit, limit), oob = scales::squish, name = "mean\nlog2FC") +
    labs(subtitle = annotation_name, x = NULL, y = NULL) +
    theme_pub(21) +
    theme(panel.grid = element_blank(),
          plot.subtitle = element_text(face = "bold", color = "black", size = rel(1.1)),
          axis.text.y = element_text(size = rel(0.85)))
}

# Assembles the (up to 4) per-category panels for one gene set into one
# figure, heights proportional to term count (same reasoning as step 13).
build_logfc_heatmap_figure <- function(ora, geneset_name, title, subtitle, limit) {
  df <- ora %>% filter(GeneSet == geneset_name) %>% add_term_logfc()
  if (nrow(df) == 0) return(invisible(NULL))

  panels <- lapply(annotation_order, function(a) make_logfc_panel(df, a, limit))
  # Wrapped LINES per panel, not raw term count - see step 13's identical
  # reasoning (make_logfc_panel wraps at the same width = 34).
  term_counts <- vapply(annotation_order, function(a) {
    terms <- unique(df$Description[df$Annotation == a])
    if (length(terms) == 0) return(0L)
    sum(vapply(terms, function(s) length(strwrap(s, width = 34)), integer(1)))
  }, integer(1))
  keep <- !vapply(panels, is.null, logical(1))
  if (!any(keep)) return(invisible(NULL))
  panels <- panels[keep]
  heights <- term_counts[keep] + 2.8

  p <- wrap_plots(panels, ncol = 1, heights = heights) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title, subtitle = subtitle, tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 25),
                    plot.subtitle = element_text(color = "grey40", size = 17),
                    plot.tag = element_text(face = "bold", size = 21))
    )
  # Figure height must scale with total row count, not be fixed - a gene
  # set with many terms across BP/CC/MF/Pathway (e.g. 20+ rows) squeezed
  # into a fixed-height image is exactly how wrapped 2-line term labels
  # end up overlapping each other. Same per-unit scaling as step 13's
  # heatmap_*_ora.png sizing (calibrated so ~50 content-units = 11in).
  list(plot = p, height = max(5, 3.0 + 0.42 * sum(heights)))
}

# One output file per gene set, matching the granularity of the
# earlier exploratory heatmaps (GO_RNA_up.png, PW_RNA_down.png, ...).
run_logfc_heatmap_set <- function(ora_file, geneset_levels, source_label, basename_prefix) {
  if (!file.exists(ora_file)) {
    message("[14_enrichment_logfc_heatmaps] Missing (skipped): ", ora_file)
    return(invisible(NULL))
  }
  ora <- read.delim(ora_file)
  ora <- add_term_logfc(ora)
  limit <- min(max(abs(c(ora$logFC_rna, ora$logFC_ribo)), na.rm = TRUE), 2)

  geneset_suffixes <- c("RNAup", "RNAdown", "RIBOup", "RIBOdown")
  for (i in seq_along(geneset_levels)) {
    built <- build_logfc_heatmap_figure(
      ora, geneset_levels[i],
      title = paste0(source_label, ": ", geneset_levels[i]),
      subtitle = "Mean log2FC of each significant term's member genes",
      limit = limit
    )
    if (is.null(built)) {
      message("[14_enrichment_logfc_heatmaps] No terms for: ", geneset_levels[i], " - skipping.")
      next
    }
    save_enrichment_fig(built$plot, paste0(basename_prefix, "_", geneset_suffixes[i]), width = 8.5, height = built$height)
  }
}

de_geneset_levels <- c("RNA-seq up-regulated", "RNA-seq down-regulated",
                        "Ribo-seq up-regulated", "Ribo-seq down-regulated")
run_logfc_heatmap_set(go_pathway_file, de_geneset_levels, "DE-based ORA", "logfc_heatmap_go_pathway")

cat("[14_enrichment_logfc_heatmaps] Done.\n")
