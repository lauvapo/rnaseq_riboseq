# ============================================================
# Step 5c: Transcription Network Enrichment ORA for DE gene sets
#
# Over-representation analysis using the transcription network
# (transcrnetwork_bustad.csv) where each regulator's target genes
# form a "regulon-like" set tested for enrichment in the
# DESeq2-significant gene sets (RNA-seq and Ribo-seq).
#
# Outputs:
#   enrichment_results/network_enrichment/network_ora_results.tsv
#   enrichment_results/network_enrichment/*.tsv (one per gene set)
# ============================================================

source(file.path(pipeline_dir, "04_enrichment_common.R"))
# (defines `all_info`, `background`, and the annotation/ORA helpers)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Load transcription network (Modulome format)
transcrnetwork <- read.csv(transcrnetwork_file, stringsAsFactors = FALSE)
cat("[05c_network_enrichment] Loaded transcription network with", nrow(transcrnetwork),
    "regulator-target edges\n")

# Create TERM2GENE table from network where:
# - TERM = Regulator (each regulator is a "term")
# - Gene = Target (targets of that regulator are the "genes")
network_term2gene <- transcrnetwork %>%
  dplyr::select(TERM = regulator_id, Gene = target_id) %>%
  distinct() %>%
  filter(!is.na(TERM), !is.na(Gene))

# Filter to genes in background (tested genes)
network_term2gene <- network_term2gene %>%
  filter(Gene %in% background)

cat("[05c_network_enrichment] Created TERM2GENE with",
    length(unique(network_term2gene$TERM)), "regulators and",
    nrow(network_term2gene), "regulator-target associations\n")

# Gene sets: same as in 05_enrichment_go_pathway.R
gene_sets <- list(
  "RNA-seq up-regulated"    = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  >  logfc_cutoff) %>% pull(name),
  "RNA-seq down-regulated"  = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  < -logfc_cutoff) %>% pull(name),
  "Ribo-seq up-regulated"   = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo >  logfc_cutoff) %>% pull(name),
  "Ribo-seq down-regulated" = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo < -logfc_cutoff) %>% pull(name)
)

# --- Run ORA for each gene set against network regulators ---
out_dir <- file.path(dir_enrichment, "network_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- lapply(names(gene_sets), function(geneset_name) {
  sig_genes <- gene_sets[[geneset_name]]
  cat("[05c_network_enrichment] ORA for:", geneset_name, "(n =", length(sig_genes), "genes)\n")

  if (length(sig_genes) == 0) {
    cat("[05c_network_enrichment]   Empty gene set, skipping\n")
    return(NULL)
  }

  ora <- enricher(
    gene = sig_genes,
    TERM2GENE = network_term2gene,
    universe = background,
    pAdjustMethod = "fdr",
    minGSSize = ora_min_gs_size,
    pvalueCutoff = ora_pvalue_cutoff,
    qvalueCutoff = ora_qvalue_cutoff
  )

  if (is.null(ora) || nrow(as.data.frame(ora)) == 0) {
    cat("[05c_network_enrichment]   No enriched regulators\n")
    return(NULL)
  }

  res_df <- as.data.frame(ora) %>%
    mutate(GeneSet = geneset_name, Annotation = "Network")

  # Save TSV
  safe_name <- gsub("[^a-zA-Z0-9_]+", "_", geneset_name)
  write.table(res_df, file.path(out_dir, paste0(safe_name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Save dotplot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- tryCatch(
      clusterProfiler::dotplot(ora, showCategory = 20) +
        ggplot2::ggtitle(paste(geneset_name, "- Network Regulators")) +
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

network_enrichment_results <- bind_rows(all_results)

# Save combined results
write.table(network_enrichment_results, file.path(dir_enrichment, "network_enrichment_ora_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[05c_network_enrichment] Done. ", nrow(network_enrichment_results),
    "enriched network regulators written to", out_dir, "\n")
