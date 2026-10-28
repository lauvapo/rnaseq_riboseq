# ============================================================
# Step 5d: iModulon ORA for DE gene sets
#
# Over-representation analysis using iModulon definitions from
# imodulons.tsv (manually curated iModulons representing
# functional modules in M. tuberculosis) for the DESeq2-significant
# gene sets (RNA-seq and Ribo-seq).
#
# Outputs:
#   enrichment_results/imodulons/imodulons_ora_results.tsv
#   enrichment_results/imodulons/*.tsv (one per gene set)
# ============================================================

source(file.path(pipeline_dir, "04_enrichment_common.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Load iModulon definitions
imodulons_file <- file.path(dir_pathways, "imodulons.tsv")
imodulons_df <- read.delim(imodulons_file, header = TRUE, stringsAsFactors = FALSE)

cat("[05d_imodulons] Loaded iModulons from", imodulons_file, "\n")
cat("[05d_imodulons] Total iModulons:", n_distinct(imodulons_df$iModulon), "\n")
cat("[05d_imodulons] Total mapped genes:", n_distinct(imodulons_df$Gene), "\n")
cat("[05d_imodulons] Total mappings:", nrow(imodulons_df), "\n")

# Convert to TERM2GENE format
# TERM = iModulon_Name (each iModulon is a "term")
# Gene = Gene (genes in that iModulon)
imodulons_term2gene <- imodulons_df %>%
  dplyr::select(TERM = iModulon_Name, Gene) %>%
  distinct() %>%
  filter(!is.na(TERM), !is.na(Gene))

# Filter to genes in background (tested genes)
imodulons_term2gene <- imodulons_term2gene %>%
  filter(Gene %in% background)

cat("[05d_imodulons] Created TERM2GENE with",
    length(unique(imodulons_term2gene$TERM)), "iModulons and",
    nrow(imodulons_term2gene), "iModulon-gene associations\n")

# Check for unmapped genes
unmapped_genes <- setdiff(imodulons_df$Gene, background)
if (length(unmapped_genes) > 0) {
  cat("[05d_imodulons] Warning:", length(unmapped_genes),
      "iModulon genes not in tested background\n")
}

# Gene sets: same as in other enrichment analyses
gene_sets <- list(
  "RNA-seq up-regulated"    = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  >  logfc_cutoff) %>% pull(name),
  "RNA-seq down-regulated"  = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  < -logfc_cutoff) %>% pull(name),
  "Ribo-seq up-regulated"   = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo >  logfc_cutoff) %>% pull(name),
  "Ribo-seq down-regulated" = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo < -logfc_cutoff) %>% pull(name)
)

# --- Run ORA for each gene set against iModulon definitions ---
out_dir <- file.path(dir_enrichment, "imodulons")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- lapply(names(gene_sets), function(geneset_name) {
  sig_genes <- gene_sets[[geneset_name]]
  cat("[05d_imodulons] ORA for:", geneset_name, "(n =", length(sig_genes), "genes)\n")

  if (length(sig_genes) == 0) {
    cat("[05d_imodulons]   Empty gene set, skipping\n")
    return(NULL)
  }

  ora <- enricher(
    gene = sig_genes,
    TERM2GENE = imodulons_term2gene,
    universe = background,
    pAdjustMethod = "fdr",
    minGSSize = ora_min_gs_size,
    pvalueCutoff = ora_pvalue_cutoff,
    qvalueCutoff = ora_qvalue_cutoff
  )

  if (is.null(ora) || nrow(as.data.frame(ora)) == 0) {
    cat("[05d_imodulons]   No enriched iModulons\n")
    return(NULL)
  }

  res_df <- as.data.frame(ora) %>%
    mutate(GeneSet = geneset_name, Annotation = "iModulon")

  # Save TSV
  safe_name <- gsub("[^a-zA-Z0-9_]+", "_", geneset_name)
  write.table(res_df, file.path(out_dir, paste0(safe_name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Save dotplot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- tryCatch(
      clusterProfiler::dotplot(ora, showCategory = 20) +
        ggplot2::ggtitle(paste(geneset_name, "- iModulon")) +
        ggplot2::theme_minimal(),
      error = function(e) NULL
    )
    if (!is.null(p)) {
      ggplot2::ggsave(file.path(out_dir, paste0(safe_name, "_dotplot.png")),
                       plot = p, width = 10, height = 6)
    }
  }

  res_df
})

imodulons_enrichment_results <- bind_rows(all_results)

# Save combined results
write.table(imodulons_enrichment_results, file.path(dir_enrichment, "imodulons_ora_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[05d_imodulons] Done. ", nrow(imodulons_enrichment_results),
    "enriched iModulons written to", out_dir, "\n")
