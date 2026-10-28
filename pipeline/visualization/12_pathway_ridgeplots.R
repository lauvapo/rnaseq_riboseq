# ============================================================
# Step 12: Ridge plots of top enriched pathway categories
#
# One ridge plot per assay (RNA-seq, Ribo-seq): a fixed, hand-picked
# set of custom-pathway ORA hits (Mycobrowser "Pathway" annotation only
# - not GO Biological Process/Cellular Component/Molecular Function;
# see SELECTED_PATHWAYS below), each ridge showing the density of
# log2FC (that assay's own logFC column) across the pathway's member
# genes. Categories are ORDERED by mean log2FC (most up-regulated at
# top), and each ridge is COLORED by -log10(p.adjust) - i.e. how
# significant that pathway's own ORA hit is - on a sequential pastel-
# to-deep purple scale (#F0E6FF light -> #6B21A8 dark).
#
# Reuses the pathway TERM2GENE table and ORA background/all_info
# defined by 04_enrichment_common.R (sourced via 05_enrichment_go_
# pathway.R, which always runs before this step in run_pipeline.R) and
# theme_pub()/save_pub() from 07_summary_figures.R - but reruns
# enricher() itself with relaxed pvalueCutoff/qvalueCutoff = 1 (instead
# of the pipeline's real 0.1/0.1 from 00_config.R), so both assays
# return a full ranked list to plot even when nothing clears the
# pipeline's formal FDR threshold for this annotation type (true for
# Ribo-seq at the time this step was written - see the "0/N pass"
# count in its own subtitle). Terms that DO clear the real pipeline
# threshold are marked with "*" on the y-axis label, so the plot never
# quietly implies more significance than 05_enrichment_go_pathway.R
# actually found.
#
# Outputs (results/):
#   ridgeplot_rna_top_pathway_categories.{png,pdf,svg}
#   ridgeplot_ribo_top_pathway_categories.{png,pdf,svg}
#   ridgeplot_combined_top_pathway_categories.{png,pdf,svg} - RNA-seq
#     and Ribo-seq as two facets sharing one pathway ordering and one
#     gene membership per pathway, for direct transcript-vs-translation
#     comparison (see that section below for details)
# ============================================================

suppressPackageStartupMessages({
  library(ggridges)
})

pathway_term2gene <- pathways2gene %>% rename(TERM = term, Gene = gene)

# --- Up/down gene sets per assay, same definition as 05_enrichment_go_pathway.R ---
gene_sets_ridge <- list(
  "RNA-seq up-regulated"    = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  >  logfc_cutoff) %>% pull(name),
  "RNA-seq down-regulated"  = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  < -logfc_cutoff) %>% pull(name),
  "Ribo-seq up-regulated"   = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo >  logfc_cutoff) %>% pull(name),
  "Ribo-seq down-regulated" = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo < -logfc_cutoff) %>% pull(name)
)

# Pathway-only ORA, relaxed thresholds so every gene set returns its
# full ranked list of tested pathways (not just the ones clearing the
# real 0.1/0.1 cutoff) - `sig_at_pipeline_cutoff` records which rows
# WOULD have passed the real cutoff, for the "*" annotation below.
run_pathway_ora_relaxed <- function(sig_genes, geneset_name) {
  if (length(sig_genes) == 0) return(NULL)
  ora <- enricher(
    gene = sig_genes,
    TERM2GENE = pathway_term2gene,
    universe = background,
    pAdjustMethod = "fdr",
    minGSSize = ora_min_gs_size,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  if (is.null(ora) || nrow(as.data.frame(ora)) == 0) return(NULL)
  as.data.frame(ora) %>%
    mutate(
      GeneSet = geneset_name,
      sig_at_pipeline_cutoff = pvalue < ora_pvalue_cutoff & qvalue < ora_qvalue_cutoff
    )
}

pathway_ora_ridge <- bind_rows(lapply(names(gene_sets_ridge), function(gs) run_pathway_ora_relaxed(gene_sets_ridge[[gs]], gs)))

# Mycobrowser pathway names carry literal HTML markup for italics
# (species/compound names, e.g. "<i>myo</i>-inositol biosynthesis") and
# a handful of HTML entities for Greek letters (e.g. "oleate
# &beta;-oxidation") - not a general HTML decoder, just what's actually
# in this annotation file. Used both to match the hand-picked pathway
# list below (which is given in already-clean form) and to label plots.
clean_pathway_name <- function(x) {
  x <- gsub("<[^>]+>", "", x)
  x <- gsub("&alpha;", "α", x, fixed = TRUE)
  x <- gsub("&beta;",  "β", x, fixed = TRUE)
  x <- gsub("&gamma;", "γ", x, fixed = TRUE)
  x <- gsub("&delta;", "δ", x, fixed = TRUE)
  x
}

# Fixed set of pathway categories to plot (all three figures below),
# rather than an automatic top-N-by-significance selection.
SELECTED_PATHWAYS <- c(
  "chorismate biosynthesis from 3-dehydroquinate",
  "mycothiol biosynthesis",
  "superpathway of L-tryptophan biosynthesis",
  "L-histidine biosynthesis",
  "mycolate biosynthesis",
  "N10-formyl-tetrahydrofolate biosynthesis",
  "mycolyl-arabinogalactan-peptidoglycan complex biosynthesis",
  "superoxide radicals degradation"
)

build_ridge_data <- function(assay_prefix, logfc_col, selected_pathways) {
  # Best (lowest) p.adjust across this assay's own up/down gene sets,
  # for the "*" marking only - independent of gene membership below.
  assay_sig <- pathway_ora_ridge %>%
    filter(grepl(paste0("^", assay_prefix), GeneSet)) %>%
    mutate(Description_clean = clean_pathway_name(Description)) %>%
    filter(Description_clean %in% selected_pathways) %>%
    group_by(Description_clean) %>%
    summarise(p.adjust = min(p.adjust), sig_at_pipeline_cutoff = any(sig_at_pipeline_cutoff), .groups = "drop")

  # Every selected pathway's FULL annotated gene membership (from
  # pathways2gene), not just the smaller "overlap with this assay's
  # significant gene set" subset an ORA hit reports - several of these
  # hand-picked pathways aren't independently significant in every
  # assay (e.g. mycothiol biosynthesis has ~0 overlap with RNA-seq's
  # significant genes), which would otherwise leave nothing to draw a
  # density curve from. Using the full membership also keeps this plot
  # and the combined one below built from identical gene sets.
  gene_rows <- tibble::tibble(Description_clean = selected_pathways) %>%
    rowwise() %>%
    reframe(
      Description_clean = Description_clean,
      gene = pathways2gene$gene[clean_pathway_name(pathways2gene$term) == Description_clean]
    )

  gene_rows %>%
    left_join(all_info[, c("name", "logFC_rna", "logFC_ribo")], by = c("gene" = "name")) %>%
    mutate(logFC = .data[[logfc_col]]) %>%
    filter(!is.na(logFC)) %>%
    left_join(assay_sig, by = "Description_clean") %>%
    mutate(
      sig_at_pipeline_cutoff = ifelse(is.na(sig_at_pipeline_cutoff), FALSE, sig_at_pipeline_cutoff),
      label_wrapped = vapply(Description_clean, function(s) paste(strwrap(s, width = 45), collapse = "\n"), character(1)),
      label_wrapped = ifelse(sig_at_pipeline_cutoff, paste0(label_wrapped, " *"), label_wrapped)
    )
}

make_ridgeplot <- function(ridge_df, assay_name) {
  # Categories ordered by log2FC (each term's mean across its member
  # genes) - geom_density_ridges draws factor level 1 at the bottom, so
  # ascending mean logFC puts the most-up-regulated term at the visual
  # top. Color is separate from order: each ridge is filled by that
  # term's own ORA -log10(p.adjust) (one solid color per ridge, not a
  # per-point gradient - geom_density_ridges_gradient silently fails to
  # render ridges backed by very few genes, several of these pathways
  # have only 2-3 member genes).
  term_means <- ridge_df %>%
    group_by(label_wrapped) %>%
    summarise(mean_logFC = mean(logFC), p.adjust = dplyr::first(p.adjust), .groups = "drop") %>%
    mutate(neglog10p = -log10(pmax(p.adjust, 1e-300)))
  ridge_df <- ridge_df %>% left_join(term_means[, c("label_wrapped", "neglog10p")], by = "label_wrapped")

  term_order <- term_means %>% arrange(mean_logFC) %>% pull(label_wrapped)
  ridge_df$label_wrapped <- factor(ridge_df$label_wrapped, levels = term_order)

  n_terms <- length(term_order)
  n_sig <- sum(grepl("\\*$", term_order))

  ggplot(ridge_df, aes(x = logFC, y = label_wrapped, fill = neglog10p)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    # rel_min_height deliberately omitted (default 0): it trims each
    # ridge's tail relative to the WHOLE PANEL's max density, not its
    # own - pathways backed by only 2-3 genes have a much lower
    # absolute peak than pathways with a dozen tightly-clustered genes,
    # so any nonzero rel_min_height here silently erases their entire
    # ridge instead of just its tail.
    geom_density_ridges(scale = 1.4, color = "white", linewidth = 0.4, alpha = 0.9,
                         jittered_points = TRUE,
                         point_shape = "|", point_size = 2.5, point_alpha = 0.6,
                         point_color = "grey20",
                         position = position_points_jitter(height = 0)) +
    scale_fill_gradient(low = "#F0E6FF", high = "#6B21A8",
                         name = expression(-log[10](p.adj))) +
    labs(
      title = paste0(assay_name, ": ", n_terms, " pathway categories"),
      subtitle = paste0("Ordered by mean log2FC, colored by -log10(p.adj) (", n_sig, "/", n_terms, " pass q<", ora_qvalue_cutoff, ", marked *)"),
      x = expression(log[2]~fold~change), y = NULL
    ) +
    theme_pub(15) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = rel(0.8))
    )
}

rna_ridge_data  <- build_ridge_data("RNA-seq",  "logFC_rna",  SELECTED_PATHWAYS)
ribo_ridge_data <- build_ridge_data("Ribo-seq", "logFC_ribo", SELECTED_PATHWAYS)

if (nrow(rna_ridge_data) > 0) {
  p_rna_ridge <- make_ridgeplot(rna_ridge_data, "RNA-seq")
  save_pub(p_rna_ridge, "ridgeplot_rna_top_pathway_categories", width = 9.5, height = 7.5)
} else {
  message("[16_pathway_ridgeplots] No RNA-seq pathway terms tested - skipping.")
}

if (nrow(ribo_ridge_data) > 0) {
  p_ribo_ridge <- make_ridgeplot(ribo_ridge_data, "Ribo-seq")
  save_pub(p_ribo_ridge, "ridgeplot_ribo_top_pathway_categories", width = 9.5, height = 7.5)
} else {
  message("[16_pathway_ridgeplots] No Ribo-seq pathway terms tested - skipping.")
}

# ============================================================
# Combined figure: RNA-seq and Ribo-seq as two panels sharing one
# pathway ordering, so a pathway's transcriptional vs. translational
# behaviour sit side by side instead of in two separate figures.
#
# Same SELECTED_PATHWAYS list as the two single-assay plots above, but
# using each pathway's FULL annotated gene membership from
# pathways2gene (not just the smaller "overlap with one significant
# gene set" subset each single-assay ORA hit reports), so the RNA-seq
# and Ribo-seq ridges for the same pathway are built from the identical
# gene list and differ only in which logFC column is plotted - a fair
# like-for-like comparison rather than two different gene subsets.
# ============================================================

# Per-(pathway, assay) significance, for the "*" annotation - best
# p.adjust across that assay's own up/down gene sets (RNA-seq and
# Ribo-seq are scored independently, unlike the shared gene list above).
assay_sig_by_pathway <- pathway_ora_ridge %>%
  mutate(Description_clean = clean_pathway_name(Description)) %>%
  filter(Description_clean %in% SELECTED_PATHWAYS) %>%
  mutate(assay = ifelse(grepl("^RNA-seq", GeneSet), "RNA-seq", "Ribo-seq")) %>%
  group_by(Description_clean, assay) %>%
  summarise(p.adjust = min(p.adjust), sig_at_pipeline_cutoff = any(sig_at_pipeline_cutoff), .groups = "drop")

combined_ridge_data <- tibble(Description_clean = SELECTED_PATHWAYS) %>%
  rowwise() %>%
  reframe(
    Description_clean = Description_clean,
    gene = pathways2gene$gene[clean_pathway_name(pathways2gene$term) == Description_clean]
  ) %>%
  left_join(all_info[, c("name", "logFC_rna", "logFC_ribo")], by = c("gene" = "name")) %>%
  tidyr::pivot_longer(c(logFC_rna, logFC_ribo), names_to = "assay", values_to = "logFC") %>%
  mutate(assay = ifelse(assay == "logFC_rna", "RNA-seq", "Ribo-seq")) %>%
  filter(!is.na(logFC)) %>%
  left_join(assay_sig_by_pathway, by = c("Description_clean", "assay")) %>%
  mutate(
    neglog10p = -log10(pmax(p.adjust, 1e-300)),
    label_wrapped = vapply(Description_clean, function(s) paste(strwrap(s, width = 40), collapse = "\n"), character(1))
  )

if (nrow(combined_ridge_data) > 0) {
  # One consistent y-order for BOTH panels: mean log2FC pooled across
  # RNA-seq and Ribo-seq together (not per-panel), so the same pathway
  # sits at the same height in both facets and is directly comparable.
  # Significance mark: "*" if that pathway passes the real pipeline
  # cutoff in EITHER assay.
  term_order_combined <- combined_ridge_data %>%
    group_by(label_wrapped) %>%
    summarise(mean_logFC_all = mean(logFC), any_sig = any(sig_at_pipeline_cutoff, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_logFC_all)

  combined_ridge_data <- combined_ridge_data %>%
    left_join(term_order_combined[, c("label_wrapped", "any_sig")], by = "label_wrapped") %>%
    mutate(label_wrapped_starred = ifelse(any_sig, paste0(label_wrapped, " *"), label_wrapped))

  term_order_starred <- term_order_combined %>%
    mutate(label_wrapped_starred = ifelse(any_sig, paste0(label_wrapped, " *"), label_wrapped)) %>%
    pull(label_wrapped_starred)

  combined_ridge_data <- combined_ridge_data %>%
    mutate(
      label_wrapped = factor(label_wrapped_starred, levels = term_order_starred),
      assay = factor(assay, levels = c("RNA-seq", "Ribo-seq"))
    )

  n_terms_combined <- length(term_order_starred)
  n_sig_combined <- sum(grepl("\\*$", term_order_starred))

  p_combined_ridge <- ggplot(combined_ridge_data, aes(x = logFC, y = label_wrapped, fill = neglog10p)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_density_ridges(scale = 1.4, color = "white", linewidth = 0.4, alpha = 0.9,
                         jittered_points = TRUE,
                         point_shape = "|", point_size = 2.2, point_alpha = 0.6,
                         point_color = "grey20",
                         position = position_points_jitter(height = 0)) +
    scale_fill_gradient(low = "#F0E6FF", high = "#6B21A8",
                         name = expression(-log[10](p.adj))) +
    facet_wrap(~assay, nrow = 1) +
    labs(
      title = paste0("RNA-seq vs Ribo-seq: ", n_terms_combined, " pathway categories"),
      subtitle = paste0("Same pathways, same gene membership, colored by each assay's own significance (", n_sig_combined, "/", n_terms_combined, " pass q<", ora_qvalue_cutoff, " in either assay, marked *)"),
      x = expression(log[2]~fold~change), y = NULL
    ) +
    theme_pub(15) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = rel(0.75)),
      strip.text = element_text(face = "bold", size = rel(1.1))
    )

  save_pub(p_combined_ridge, "ridgeplot_combined_top_pathway_categories", width = 13, height = 7.5)
} else {
  message("[16_pathway_ridgeplots] No pathway terms available for the combined plot - skipping.")
}

cat("[16_pathway_ridgeplots] Done.\n")
