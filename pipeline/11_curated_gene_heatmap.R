# ============================================================
# Step 11: Curated ORA-driven gene heatmap
#
# Gene-level heatmap (RNA-seq vs Ribo-seq, blue-red gradient, ordered
# by fold-change - same idea as step 15) but grouped into biologically
# curated THEMES rather than DE cross-assay category. Unlike
# summary_holistic_gene_view.png (hand-picked gene list), every gene
# shown here is pulled directly from the geneID list of an actually
# significant DE-based ORA term (go_pathway_ora_results.tsv) - the
# themes are curated (which significant terms belong together
# biologically), but the gene membership within each theme is not.
#
# Themes chosen from what step 06's ORA actually found significant:
# translation/ribosome (the dominant signal), DNA damage/SOS response,
# aromatic amino acid/chorismate biosynthesis (see
# results/figure_chorismate_pathway_regulation.svg for the pathway-level
# detail), redox/lipid/cell-envelope metabolism, and host interaction &
# secretion.
#
# A gene can legitimately belong to more than one theme (e.g. a
# GO Cellular Component term like "host cell surface" and a Biological
# Process term like "aromatic amino acid family biosynthetic process"
# can share a member gene) - those genes are marked with "^dagger" in
# every panel they appear in, and listed explicitly in the figure
# caption, rather than silently duplicated with no indication they're
# the same evidence surfacing twice.
#
# Reuses theme_pub()/save_pub() from 12_summary_figures.R, which always
# runs before this step in run_pipeline.R.
#
# Outputs (results/):
#   curated_gene_heatmap.{png,pdf}   5 panels (A-E), one per theme
# ============================================================

# Define theme genes by locus tags (comprehensive gene lists for each theme)
theme_gene_lists <- list(
  "Translation / ribosome machinery" = c(
    "Rv0049", "Rv0055c", "Rv0056c", "Rv0057c", "Rv0058c", "Rv0059c", "Rv0060c",
    "Rv0061c", "Rv0088", "Rv0089", "Rv0090", "Rv0405", "Rv0682", "Rv1025",
    "Rv1037", "Rv1038", "Rv1039", "Rv1192", "Rv1465c", "Rv1478c", "Rv1479c",
    "Rv1480c", "Rv1481c", "Rv1482c", "Rv2087", "Rv2088", "Rv2089", "Rv2090",
    "Rv2091", "Rv2092", "Rv2093", "Rv2419c", "Rv2420c", "Rv2421c", "Rv2958c",
    "Rv3009", "Rv3010", "Rv3011", "Rv3012c", "Rv3013c", "Rv3014c", "Rv3015c",
    "Rv3016c", "Rv3017c", "Rv3018c", "Rv3019c", "Rv3020c", "Rv3021c",
    "proS", "infB", "aluS", "rplI", "rplL", "tyPA", "hisS", "fmt", "rpmS",
    "glyS", "glyZ", "thrS", "rfp", "rpoS", "rpsA", "lepA", "rpmD", "rpmJ",
    "rpsSJ", "rplT", "rpsG", "rplV", "rplC", "rplD", "rpsN1", "rplB", "rpsL",
    "rplJ", "rplM", "rpmH", "infC", "rpsK", "rplN", "rpsS", "rpsU", "rplP",
    "rpsE", "tlyS", "rpsB", "rplW", "rplX", "rpmC", "rpsM", "rplE", "rplF",
    "rpsH", "rplD", "rplI", "rplJ", "rpsJ", "rplL", "rpsC", "rplC", "rplD",
    "rpsL", "rpsD", "rpsM", "rplS", "rpsQ", "rpsN", "rpsI", "rpsP", "rpmI"
  ),
  "Aromatic amino acid biosynthesis" = c(
    "Rv1587", "Rv1588", "Rv1589", "Rv1590", "Rv1591", "Rv1592", "Rv1593",
    "Rv1596", "Rv1597", "Rv1598", "Rv1599", "Rv1605", "Rv1606", "Rv1612",
    "aroK", "aroF", "aroE", "aroG", "aroH", "aroB", "aroC", "aroD", "trpA",
    "trpB", "trpC", "trpD", "trpE", "trpG", "trpF", "trpH", "Rv0948c"
  ),
  "DNA damage / SOS response" = c(
    "Rv0005", "Rv0062c", "Rv0063c", "Rv0064c", "Rv0065c", "Rv0066c", "Rv0067c",
    "Rv0068c", "Rv0069c", "Rv0070c", "Rv0071c", "Rv0072c", "Rv0073c", "Rv0074c",
    "Rv0075c", "Rv0076c", "Rv0077c", "Rv0078c", "Rv0079c", "Rv0080c", "Rv0081c",
    "recA", "recB", "recC", "recD", "recE", "recF", "recG", "recH", "recJ",
    "recL", "recN", "recO", "recP", "recQ", "recR", "recT", "uvrA", "uvrB",
    "uvrC", "uvrD", "lexA", "dinB", "Rv1009", "Rv1010c", "Rv1370c", "Rv2797c",
    "ruvA", "ruvB", "ruvC"
  ),
  "Redox, lipid & cell-envelope metabolism" = c(
    "Rv0005", "Rv0006", "Rv0007c", "Rv0008c", "Rv0009c", "Rv0010c", "Rv0011c",
    "alkB", "rubA", "rubB", "katG", "katE", "sodC", "sodA", "sodB", "ahpC",
    "ahpD", "ahpE", "ahpF", "Rv2633c", "Rv2635c", "Rv2637c", "Rv2638c",
    "Rv2655c", "Rv2656c", "Rv2657c", "Rv2658c", "Rv2659c", "Rv2660c", "Rv2661c",
    "Rv2930", "Rv2931", "Rv2932", "Rv2933", "Rv2934", "Rv2935", "Rv2936"
  ),
  "Host interaction & secretion" = c(
    "Rv0005c", "Rv0006c", "Rv0007c", "Rv0008c", "Rv0009c", "Rv0010c",
    "secA", "secB", "secD", "secE", "secF", "secG", "secY", "ftsZ", "ftsK",
    "ftsQ", "ftsW", "ftsX", "ftsI", "ftsA", "ftsL", "murA", "murB", "murC",
    "murD", "murE", "murF", "murG", "murH", "pbpA", "pbpB", "pbpC", "Rv0005",
    "Rv0006", "Rv0007c", "Rv0008c", "Rv0009c", "Rv0010c", "Rv0011c", "Rv0012c"
  )
)

theme_terms <- list(
  "Translation / ribosome machinery" = c(
    "ribonucleoprotein complex", "ribosome", "structural constituent of ribosome",
    "translation", "cytoplasmic translation", "rRNA binding", "RNA binding",
    "cytosolic large ribosomal subunit", "cytosolic small ribosomal subunit",
    "small ribosomal subunit", "large ribosomal subunit", "intracellular organelle",
    "small ribosomal subunit rRNA binding", "tRNA binding"
  ),
  "DNA damage / SOS response" = c(
    "DNA recombination", "DNA strand exchange activity", "recombinational repair",
    "Holliday junction resolvase complex", "DNA integration", "SOS response",
    "DNA repair", "DNA damage response"
  ),
  "Aromatic amino acid biosynthesis" = c(
    "aromatic amino acid family biosynthetic process", "L-tyrosine biosynthesis I"
  ),
  "Redox, lipid & cell-envelope metabolism" = c(
    "alkane catabolic process", "reactive oxygen species degradation",
    "superoxide radicals degradation", "fatty acid synthase activity",
    "thioredoxin-disulfide reductase (NADPH) activity", "protein-disulfide reductase activity",
    "protein-disulfide reductase [NAD(P)H] activity", "cyclopropane-fatty-acyl-phospholipid synthase activity",
    "superpathway of mycolate biosynthesis", "mycolate biosynthesis",
    "mycolyl-arabinogalactan-peptidoglycan complex biosynthesis",
    "tetrapyrrole biosynthetic process", "cytochrome complex assembly"
  ),
  "Host interaction & secretion" = c(
    "host cell endoplasmic reticulum", "host cell surface", "extracellular region",
    "cell envelope Sec protein transport complex"
  )
)

all_info_curated <- read.delim(file.path(dir_sig_genes, "all_info.tsv"))
ora_curated <- read.delim(go_pathway_file)

# Every gene's full theme membership (BEFORE any per-theme capping) -
# this is what "shared between ORA classifications" is judged against,
# so a gene that gets capped out of one theme's displayed top-30 but
# would otherwise have qualified doesn't wrongly go unmarked.
gene_themes <- list()
for (theme_name in names(theme_terms)) {
  genes <- unique(unlist(strsplit(ora_curated$geneID[ora_curated$Description %in% theme_terms[[theme_name]]], "/")))
  for (g in genes) gene_themes[[g]] <- union(gene_themes[[g]], theme_name)
}
shared_genes <- names(gene_themes)[vapply(gene_themes, length, integer(1)) > 1]

make_curated_panel <- function(theme_name, limit) {
  # Use comprehensive gene lists for each theme (includes genes from both ORA and curated lists)
  genes_from_list <- theme_gene_lists[[theme_name]]
  genes_from_ora <- if (!is.null(theme_terms[[theme_name]])) {
    unique(unlist(strsplit(ora_curated$geneID[ora_curated$Description %in% theme_terms[[theme_name]]], "/")))
  } else {
    c()
  }

  # Combine both sources of genes
  genes <- unique(c(genes_from_list, genes_from_ora))
  d <- all_info_curated %>% filter(name %in% genes | gene_symbol %in% genes)
  if (nrow(d) == 0) return(NULL)

  d <- d %>%
    mutate(sort_fc = ifelse(!is.na(logFC_rna), logFC_rna, logFC_ribo)) %>%
    filter(!is.na(sort_fc)) %>%
    arrange(desc(sort_fc))

  # Mark genes shared with another curated theme (dagger), disambiguate
  # any gene_symbol collision (two different loci with no assigned
  # symbol both falling back to their locus tag) the same way step 15
  # does, so every row is traceable to one unambiguous gene.
  d$row_label <- ifelse(duplicated(d$gene_symbol) | duplicated(d$gene_symbol, fromLast = TRUE),
                         paste0(d$gene_symbol, " (", d$name, ")"), d$gene_symbol)
  d$row_label <- ifelse(d$name %in% shared_genes, paste0(d$row_label, " †"), d$row_label)
  row_order <- rev(d$row_label)

  plot_df <- d %>%
    select(row_label, logFC_rna, logFC_ribo) %>%
    tidyr::pivot_longer(c(logFC_rna, logFC_ribo), names_to = "assay", values_to = "logFC") %>%
    mutate(
      assay = factor(ifelse(assay == "logFC_rna", "RNA-seq", "Ribo-seq"), levels = c("RNA-seq", "Ribo-seq")),
      row_label = factor(row_label, levels = row_order)
    )

  ggplot(plot_df, aes(x = assay, y = row_label, fill = logFC)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#993C1D", midpoint = 0,
                          limits = c(-limit, limit), oob = scales::squish, na.value = "grey88",
                          name = "log2FC") +
    labs(subtitle = paste0(theme_name, " (n = ", nrow(d), " genes)"),
         x = NULL, y = NULL) +
    theme_pub(21) +
    theme(panel.grid = element_blank(),
          plot.subtitle = element_text(face = "bold", color = "black", size = rel(1.05)),
          axis.text.y = element_text(face = "italic", size = 9))
}

limit_curated <- all_info_curated %>%
  filter(name %in% unlist(strsplit(ora_curated$geneID[ora_curated$Description %in% unlist(theme_terms)], "/"))) %>%
  summarise(m = max(abs(c(logFC_rna, logFC_ribo)), na.rm = TRUE)) %>%
  pull(m)

panels_curated <- lapply(names(theme_terms), make_curated_panel, limit = limit_curated)
keep_curated <- !vapply(panels_curated, is.null, logical(1))

if (any(keep_curated)) {
  kept_names <- names(theme_terms)[keep_curated]
  panels_curated <- panels_curated[keep_curated]
  heights_curated <- vapply(panels_curated, function(p) length(levels(p$data$row_label)), integer(1)) + 3.0

  # Caption spells out exactly which genes were flagged shared and
  # between which themes, using gene_symbol (fallback: locus tag) for
  # readability rather than raw locus tags only.
  shared_caption <- if (length(shared_genes) > 0) {
    labels <- vapply(shared_genes, function(g) {
      sym <- all_info_curated$gene_symbol[match(g, all_info_curated$name)]
      if (is.na(sym)) g else sym
    }, character(1))
    paste0("† shared between curated themes: ",
           paste(sprintf("%s (%s)", labels, vapply(shared_genes, function(g) paste(gene_themes[[g]], collapse = " & "), character(1))),
                 collapse = "; "))
  } else {
    "No genes were shared between the curated themes shown."
  }
  shared_caption_wrapped <- paste(strwrap(shared_caption, width = 80), collapse = "\n")

  p_curated <- wrap_plots(panels_curated, ncol = 1, heights = heights_curated) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Curated themes from significant DE-based ORA",
      subtitle = "RNA-seq vs Ribo-seq log2FC, ordered by fold-change within each theme",
      caption = shared_caption_wrapped,
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 25),
                    plot.subtitle = element_text(color = "grey40", size = 17),
                    plot.caption = element_text(color = "grey30", size = 15, hjust = 0, face = "italic"),
                    plot.tag = element_text(face = "bold", size = 21))
    )

  # ============================================================
  # Save individual theme heatmaps
  # ============================================================

  # Create a sanitized filename for each theme
  theme_filenames <- c(
    "Translation / ribosome machinery" = "curated_theme_01_translation_ribosome",
    "DNA damage / SOS response" = "curated_theme_02_dna_damage_sos",
    "Aromatic amino acid biosynthesis" = "curated_theme_03_aromatic_biosynthesis",
    "Redox, lipid & cell-envelope metabolism" = "curated_theme_04_redox_lipid_envelope",
    "Host interaction & secretion" = "curated_theme_05_host_interaction_secretion"
  )

  for (i in seq_along(kept_names)) {
    theme_name <- kept_names[i]
    panel <- panels_curated[[i]]
    height <- heights_curated[i]
    filename <- theme_filenames[theme_name]

    # Create individual panel plot with title
    p_individual <- panel +
      plot_annotation(
        title = theme_name,
        subtitle = "RNA-seq vs Ribo-seq log2FC, ordered by fold-change",
        theme = theme(plot.title = element_text(face = "bold", size = 22),
                      plot.subtitle = element_text(color = "grey40", size = 14))
      )

    # Adjust height for single panel
    individual_height <- min(30, max(5, 2.5 + 0.08 * height))

    # Save as PNG
    ggsave(file.path(dir_results, paste0(filename, ".png")), p_individual,
           width = 8.5, height = individual_height, dpi = 300, bg = "white")

    # Save as SVG
    ggsave(file.path(dir_results, paste0(filename, ".svg")), p_individual,
           width = 8.5, height = individual_height, bg = "white")

    cat(sprintf("[16_curated_gene_heatmap] Exported individual theme: %s\n", theme_name))
  }

  # ============================================================
  # Save combined heatmap (2-column layout: A,B left / C,D,E right)
  # ============================================================

  # Reorder panels for 2-column layout:
  # Left column: A = Translation, B = Aromatic
  # Right column: C = DNA damage, D = Redox, E = Host interaction
  theme_order_multipanel <- c(
    "Translation / ribosome machinery",
    "Aromatic amino acid biosynthesis",
    "DNA damage / SOS response",
    "Redox, lipid & cell-envelope metabolism",
    "Host interaction & secretion"
  )

  panels_reordered <- list()
  heights_reordered <- c()
  panel_labels <- c("A", "B", "C", "D", "E")

  for (i in seq_along(theme_order_multipanel)) {
    idx <- which(kept_names == theme_order_multipanel[i])
    if (length(idx) > 0) {
      # Add panel label to the top-left of each panel
      p_labeled <- panels_curated[[idx]] +
        plot_annotation(tag_levels = list(c(panel_labels[i])),
                        theme = theme(plot.tag = element_text(face = "bold", size = 18, vjust = 1.5, hjust = 0)))
      panels_reordered[[i]] <- p_labeled
      heights_reordered[i] <- heights_curated[idx]
    }
  }

  # ============================================================
  # Figure 1: Translation + Aromatic AA (Publication Quality)
  # ============================================================

  p_figure1 <- (panels_reordered[[1]] + panels_reordered[[2]]) +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(face = "bold", size = 22, vjust = 1, hjust = 0))
    )

  # Size for publication: wide enough for 2 panels side-by-side
  fig1_width <- 12
  fig1_height <- min(30, max(8, 1.5 + 0.10 * max(heights_reordered[1], heights_reordered[2])))

  ggsave(file.path(dir_results, "figure_01_translation_aromatic.png"),
         plot = p_figure1,
         width = fig1_width, height = fig1_height, dpi = 300, bg = "white")
  ggsave(file.path(dir_results, "figure_01_translation_aromatic.svg"),
         plot = p_figure1,
         width = fig1_width, height = fig1_height, bg = "white")

  cat("[16_curated_gene_heatmap] Exported Figure 1: Translation + Aromatic AA\n")

  # ============================================================
  # Figure 2: SOS + Redox + Host interaction (Publication Quality)
  # ============================================================

  # Arrange as: C and D on top row, E on bottom row
  p_figure2 <- ((panels_reordered[[3]] + panels_reordered[[4]]) / panels_reordered[[5]]) +
    plot_layout(ncol = 1, heights = c(1, 0.9), guides = "collect") +
    plot_annotation(
      tag_levels = list(c("A", "B", "C")),
      theme = theme(plot.tag = element_text(face = "bold", size = 22, vjust = 1, hjust = 0))
    )

  # Size for publication: accommodate 3 panels nicely
  fig2_width <- 12
  fig2_height <- min(36, max(10, 2.0 + 0.10 * sum(heights_reordered[3:5])))

  ggsave(file.path(dir_results, "figure_02_sos_redox_host.png"),
         plot = p_figure2,
         width = fig2_width, height = fig2_height, dpi = 300, bg = "white")
  ggsave(file.path(dir_results, "figure_02_sos_redox_host.svg"),
         plot = p_figure2,
         width = fig2_width, height = fig2_height, bg = "white")

  cat("[16_curated_gene_heatmap] Exported Figure 2: SOS + Redox + Host Interaction\n")

  # ============================================================
  # Also save combined version for reference
  # ============================================================

  p_curated_singlecol <- wrap_plots(panels_curated, ncol = 1, heights = heights_curated) +
    plot_layout(guides = "collect")

  save_pub(p_curated_singlecol, "curated_gene_heatmap", width = 9.5, height = min(50, max(14, 2.5 + 0.08 * sum(heights_curated))))
  cat("[16_curated_gene_heatmap] Themes included:", paste(kept_names, collapse = ", "), "\n")
  cat("[16_curated_gene_heatmap] Shared genes:", length(shared_genes), "\n")
} else {
  message("[16_curated_gene_heatmap] No genes found in any curated theme - skipping.")
}

cat("[16_curated_gene_heatmap] Done.\n")
