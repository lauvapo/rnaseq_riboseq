# ============================================================
# Step 8: Clustered ORA heatmaps
#
# summary_enrichment_overview.png (step 12) shows how many enriched
# terms landed in each gene-set x annotation-category cell, but not
# WHICH terms or how they relate to each other. This step goes one
# level deeper: for the DE-based ORA results and each annotation
# category (Biological Process, Cellular Component, Molecular
# Function, Pathway), builds a term x gene-set heatmap of
# -log10(adjusted p), with terms hierarchically clustered so terms that
# are enriched in the same gene sets end up next to each other instead
# of in arbitrary/alphabetical order.
#
# Reuses theme_pub() defined in 12_summary_figures.R, which always runs
# immediately before this step in run_pipeline.R (save_enrichment_fig()
# below is this step's own save routine, since these outputs belong in
# dir_enrichment rather than dir_results where save_pub() writes).
#
# Outputs (enrichment_results/):
#   heatmap_go_pathway_ora.{png,pdf}   DE-based ORA, 4 panels (A-D)
# ============================================================

annotation_order <- c("Biological Process", "Cellular Component", "Molecular Function", "Pathway")

# save_pub() from step 12 always writes into dir_results; these heatmaps
# belong with the rest of the ORA detail output in dir_enrichment
# instead, so this is its own small save routine rather than repurposing
# save_pub() with a relative-path hack.
save_enrichment_fig <- function(plot, basename, width, height) {
  ggsave(file.path(dir_enrichment, paste0(basename, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(dir_enrichment, paste0(basename, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(dir_enrichment, paste0(basename, ".svg")), plot,
         width = width, height = height, bg = "white")
}

# Hierarchically clusters terms (rows) by their -log10(p.adjust) profile
# across gene sets (columns), via average-linkage clustering on
# Euclidean distance - the standard way to make a term x condition
# heatmap readable instead of leaving rows in arbitrary order. With
# fewer than 3 terms there's nothing meaningful to cluster.
cluster_term_order <- function(mat) {
  if (nrow(mat) < 3) return(rownames(mat))
  rownames(mat)[hclust(dist(mat), method = "average")$order]
}

# Builds one clustered heatmap panel for a single annotation category.
# `geneset_levels` fixes the column order (and includes columns with no
# hits for this category, shown as blank/white tiles - absence of
# enrichment is as informative here as presence); `geneset_labels` is
# the shorter text shown on the panel's x-axis. `cap` is the shared
# color-scale ceiling computed once per figure (see build_ora_heatmap_
# figure) - without a shared cap, ribosome-related terms with
# -log10(p.adjust) up to ~30 make every other, still-significant term
# (typically 1-8) visually indistinguishable from zero on a linear
# scale, AND each panel would auto-scale to its own range independently
# so the four panels' colors wouldn't even mean the same thing.
make_ora_heatmap_panel <- function(df, geneset_levels, geneset_labels, annotation_name, cap) {
  df <- df %>% filter(Annotation == annotation_name)
  if (nrow(df) == 0) return(NULL)

  df <- df %>%
    mutate(GeneSet = factor(GeneSet, levels = geneset_levels)) %>%
    group_by(Description, GeneSet) %>%
    slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
    ungroup()

  wide <- df %>%
    select(Description, GeneSet, p.adjust) %>%
    tidyr::complete(Description, GeneSet = factor(geneset_levels, levels = geneset_levels)) %>%
    mutate(neglog10 = ifelse(is.na(p.adjust), 0, -log10(pmax(p.adjust, 1e-300)))) %>%
    select(-p.adjust) %>%
    tidyr::pivot_wider(names_from = GeneSet, values_from = neglog10) %>%
    as.data.frame()
  rownames(wide) <- wide$Description
  mat <- as.matrix(wide[, -1, drop = FALSE])

  term_order <- cluster_term_order(mat)

  plot_df <- df %>%
    tidyr::complete(Description, GeneSet = factor(geneset_levels, levels = geneset_levels)) %>%
    mutate(
      neglog10 = ifelse(is.na(p.adjust), 0, -log10(pmax(p.adjust, 1e-300))),
      Description = factor(Description, levels = rev(term_order)),
      label_wrapped = vapply(as.character(Description),
                              function(s) paste(strwrap(s, width = 34), collapse = "\n"), character(1)),
      GeneSet = factor(GeneSet, levels = geneset_levels, labels = geneset_labels)
    )
  # strwrap() above was computed per-row from the (unordered) character
  # vector, so label_wrapped needs its own factor built from the same
  # rev(term_order) sequence, wrapped consistently, for the y-axis.
  wrapped_levels <- vapply(rev(term_order), function(s) paste(strwrap(s, width = 34), collapse = "\n"), character(1))
  plot_df$label_wrapped <- factor(plot_df$label_wrapped, levels = wrapped_levels)

  ggplot(plot_df, aes(x = GeneSet, y = label_wrapped, fill = neglog10)) +
    geom_tile(color = "white", linewidth = 0.6) +
    scale_fill_gradient(low = "grey96", high = "#0F6E56", limits = c(0, cap),
                         oob = scales::squish, name = expression(-log[10](p.adj))) +
    labs(subtitle = annotation_name, x = NULL, y = NULL) +
    theme_pub(21) +
    theme(panel.grid = element_blank(),
          plot.subtitle = element_text(face = "bold", color = "black", size = rel(1.1)),
          axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = rel(0.85)))
}

# Assembles the 4 per-category panels (that actually have terms) into
# one figure, heights proportional to term count - a busy 15-term
# Molecular Function panel would otherwise get the same height as a
# 1-term Pathway panel and cram its labels illegibly.
build_ora_heatmap_figure <- function(ora_file, geneset_levels, geneset_labels, title, subtitle) {
  if (!file.exists(ora_file)) {
    message("[13_enrichment_heatmaps] Missing (skipped): ", ora_file)
    return(invisible(NULL))
  }
  ora <- read.delim(ora_file)

  # Shared color-scale ceiling for every panel in this figure: the
  # source's own max significance, but never above 10 (p.adjust = 1e-10
  # is already about as significant as this kind of test gets, so
  # letting a handful of extreme ribosome-related terms push the cap to
  # 30+ would wash out every more-modest-but-real signal - see the
  # comment on make_ora_heatmap_panel). A source whose own max never
  # reaches 10 is left uncapped at its true max, so it isn't stretched
  # across a mostly-empty range either.
  cap <- min(max(-log10(pmax(ora$p.adjust, 1e-300))), 10)

  panels <- lapply(annotation_order, function(a) make_ora_heatmap_panel(ora, geneset_levels, geneset_labels, a, cap))
  # Counts wrapped LINES per panel, not just terms - a panel full of
  # long GO term names that each wrap to 2 lines (width = 34, matching
  # make_ora_heatmap_panel's strwrap()) needs roughly twice the height
  # of a same-term-count panel of short names, or its rows visibly
  # collide regardless of how generous the overall figure height is.
  term_counts <- vapply(annotation_order, function(a) {
    terms <- unique(ora$Description[ora$Annotation == a])
    if (length(terms) == 0) return(0L)
    sum(vapply(terms, function(s) length(strwrap(s, width = 34)), integer(1)))
  }, integer(1))
  keep <- !vapply(panels, is.null, logical(1))
  if (!any(keep)) {
    message("[13_enrichment_heatmaps] No enriched terms in any category - skipping: ", ora_file)
    return(invisible(NULL))
  }
  panels <- panels[keep]
  heights <- term_counts[keep] + 3.0

  p <- wrap_plots(panels, ncol = 1, heights = heights) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title, subtitle = subtitle, tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 25),
                    plot.subtitle = element_text(color = "grey40", size = 17),
                    plot.tag = element_text(face = "bold", size = 21))
    )
  # Height must scale with total wrapped-line count, not be fixed - a
  # fixed total height means shifting more units into a busy panel only
  # steals inches from the other panels (or, worse, shrinks everyone's
  # per-unit share), it doesn't add real space. Same per-unit scaling as
  # steps 14-16.
  list(plot = p, height = max(5, 3.0 + 0.4 * sum(heights)))
}

de_geneset_levels <- c("RNA-seq up-regulated", "RNA-seq down-regulated",
                        "Ribo-seq up-regulated", "Ribo-seq down-regulated")
de_geneset_labels <- c("RNA-seq\nup", "RNA-seq\ndown", "Ribo-seq\nup", "Ribo-seq\ndown")

built_de <- build_ora_heatmap_figure(
  go_pathway_file, de_geneset_levels, de_geneset_labels,
  "Clustered enrichment: DE-based ORA",
  "Terms hierarchically clustered by significance pattern across gene sets"
)
if (!is.null(built_de)) {
  save_enrichment_fig(built_de$plot, "heatmap_go_pathway_ora", width = 11.5, height = built_de$height)
}

cat("[13_enrichment_heatmaps] Done.\n")
