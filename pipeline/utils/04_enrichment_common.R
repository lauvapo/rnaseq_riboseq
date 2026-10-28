# ============================================================
# Step 4: Shared GO / pathway ORA machinery
#
# ORA (over-representation analysis) asks: out of the genes annotated
# with a given GO term or pathway, is a bigger fraction of them in my
# gene-of-interest list than you'd expect by chance, given how many
# genes total have that annotation? This is a hypergeometric test
# (clusterProfiler::enricher() below), run separately for every GO
# term / pathway against every gene set (RNA-seq up, RNA-seq down,
# etc.), with the p-values FDR-corrected across all terms tested
# within each run.
#
# Builds the TERM2GENE tables used by every enrichment step (06, 07, 09)
# from the local copies of the Mycobrowser GO-term and pathway gene
# annotation exports, and defines the enrichment-and-save helper they
# share. Not meant to be run standalone; sourced by 06, 07, and 09.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(janitor)
})

all_info <- read.delim(file.path(dir_sig_genes, "all_info.tsv"))
# The "background"/"universe" for ORA: every gene that was actually
# tested (i.e. could in principle have ended up in a significant gene
# set), NOT the whole genome. Using the whole genome as background
# would overstate enrichment, since genes filtered out for low counts
# never had a chance to appear in any gene set being tested.
background <- unique(na.omit(all_info$name))

# --- GO terms (Biological Process / Cellular Component / Molecular Function) ---
# go_annotation_file lists, per gene SYMBOL, which GO terms apply in
# each of the 3 GO namespaces (one column per namespace, multiple terms
# per cell separated by " // "). left_join(gene_map) adds each gene's
# locus tag, since that's the identifier used everywhere else in this
# pipeline (symbols aren't guaranteed unique/available - see 00_config.R).
go_annot <- read.delim(go_annotation_file, header = TRUE, stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  left_join(gene_map, by = c("genes" = "name"))

# clusterProfiler::enricher() expects a TERM2GENE table: one row per
# (term, gene) pair. separate_rows() explodes the multi-term cells
# above into one row per term, so a gene annotated with 3 GO terms
# becomes 3 rows.
create_term2gene <- function(df, term_col) {
  df %>%
    dplyr::select(TERM = {{ term_col }}, Gene = locus_tag) %>%
    separate_rows(TERM, sep = " // ") %>%
    filter(TERM != "" & !is.na(Gene))
}

term2gene_list <- list(
  "Biological Process" = create_term2gene(go_annot, go_terms_biological_process),
  "Cellular Component"  = create_term2gene(go_annot, go_terms_cellular_component),
  "Molecular Function"  = create_term2gene(go_annot, go_terms_molecular_function)
)

# --- Custom pathways ---
pathway_annot <- read.delim(pathway_annotation_file, header = TRUE, stringsAsFactors = FALSE)

pathways2gene <- pathway_annot %>%
  separate_rows(Genes.of.pathway, sep = " // ") %>%
  dplyr::select(term = Pathways, gene = Genes.of.pathway) %>%
  left_join(gene_map, by = c("gene" = "name")) %>%
  mutate(gene = ifelse(is.na(locus_tag), gene, locus_tag)) %>%
  dplyr::select(term, gene) %>%
  distinct()

# Some pathway gene symbols have no entry in gene_map at all (even after
# the Mycobrowser-v5 NA patch above), so the left_join fallback leaves them
# as a raw symbol instead of a locus tag. Those rows can never match
# `background` or any sig_genes vector (both are locus tags), so they are
# silently invisible to every ORA test from here on. Not fixable here -
# gene_map is an external file this pipeline doesn't own - but surfaced
# instead of failing silently, so a genuinely relevant pathway gene missing
# from an enrichment result can be traced back to this cause rather than
# assumed to be biological. Valid locus tags are either "Rv"-prefixed
# (protein-coding genes) or "MTB0000xx" (Mycobrowser's ID for ncRNAs/tRNAs).
unmapped_pathway_genes <- pathways2gene %>%
  filter(!grepl("^(Rv[0-9]+[A-Za-z]?|MTB[0-9]+)$", gene)) %>%
  pull(gene) %>%
  unique() %>%
  sort()
if (length(unmapped_pathway_genes) > 0) {
  message("[enrichment] ", length(unmapped_pathway_genes),
          " pathway gene symbols have no locus tag in gene_map and will ",
          "never match the ORA universe/gene sets: ",
          paste(unmapped_pathway_genes, collapse = ", "))
}

# The 4 annotation types every gene set gets tested against: 3 GO
# namespaces plus the custom pathway list.
annotation_sets <- c(term2gene_list, list("Pathway" = pathways2gene %>% rename(TERM = term, Gene = gene)))

# --- Run ORA for one gene set against one annotation, save plot + table ---
# Wraps clusterProfiler::enricher() (the actual hypergeometric test)
# with bookkeeping: skip cleanly if the gene set is empty or nothing
# comes out significant (both are expected outcomes, not errors), and
# save whatever DOES come out as both a TSV and a dotplot PNG.
run_and_save_ora <- function(sig_genes, term2gene, annotation_name, geneset_name, out_dir) {
  if (length(sig_genes) == 0) {
    message("[enrichment] Empty gene set, skipping: ", geneset_name, " / ", annotation_name)
    return(NULL)
  }

  ora <- enricher(
    gene = sig_genes,
    TERM2GENE = term2gene,
    universe = background,
    pAdjustMethod = "fdr",
    minGSSize = ora_min_gs_size,
    pvalueCutoff = ora_pvalue_cutoff,
    qvalueCutoff = ora_qvalue_cutoff
  )

  if (is.null(ora) || nrow(as.data.frame(ora)) == 0) {
    message("[enrichment] No enriched terms: ", geneset_name, " / ", annotation_name)
    return(NULL)
  }

  res_df <- as.data.frame(ora) %>%
    mutate(GeneSet = geneset_name, Annotation = annotation_name)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  safe_name <- gsub("[^a-zA-Z0-9_]+", "_", paste(geneset_name, annotation_name, sep = "_"))
  write.table(res_df, file.path(out_dir, paste0(safe_name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- tryCatch(
      clusterProfiler::dotplot(ora, showCategory = 20) +
        ggplot2::ggtitle(paste(geneset_name, "-", annotation_name)) +
        ggplot2::theme_minimal(),
      error = function(e) NULL
    )
    if (!is.null(p)) {
      ggplot2::ggsave(file.path(out_dir, paste0(safe_name, "_dotplot.png")),
                       plot = p, width = 10, height = 6)
    }
  }

  res_df
}

# --- Run ORA for a gene set against every annotation (GO x3 + Pathway) ---
# The entry point steps 06/07/09 actually call: one gene set in, one
# combined results table out (rows from all 4 annotation types stacked
# together, or an empty result if nothing was significant anywhere).
run_ora_all_annotations <- function(sig_genes, geneset_name, out_dir) {
  results <- lapply(names(annotation_sets), function(annotation_name) {
    run_and_save_ora(sig_genes, annotation_sets[[annotation_name]], annotation_name, geneset_name, out_dir)
  })
  bind_rows(results)
}
