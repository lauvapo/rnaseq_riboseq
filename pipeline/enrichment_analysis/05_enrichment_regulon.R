# ============================================================
# Step 5b: Regulon ORA for DE gene sets
#
# Over-representation analysis using regulon definitions from
# selected_regulon.gmt (manually curated transcription factor
# targets) for the DESeq2-significant gene sets (RNA-seq and Ribo-seq).
#
# Outputs:
#   enrichment_results/regulon/regulon_ora_results.tsv
#   enrichment_results/regulon/*.tsv (one per gene set)
# ============================================================

source(file.path(pipeline_dir, "utils/04_enrichment_common.R"))
# (defines `all_info`, `background`, and the annotation/ORA helpers)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# --- Parse GMT file (regulon definitions) ---
# Format: Regulon_Name \t description \t Gene1 \t Gene2 \t ...
parse_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)
  regulon_list <- list()

  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      regulon_name <- parts[1]
      # Skip column 2 (description), columns 3+ are genes
      genes <- parts[3:length(parts)]
      genes <- genes[genes != ""]
      regulon_list[[regulon_name]] <- genes
    }
  }

  regulon_list
}

# Read GMT and convert to TERM2GENE format
regulons <- parse_gmt(regulon_annotation_file)
cat("[05b_regulon] Loaded", length(regulons), "regulons from GMT file\n")

# Convert to clusterProfiler-compatible TERM2GENE table
# GMT already contains locus tags, but verify they're valid and filter background
regulon_term2gene <- lapply(names(regulons), function(reg_name) {
  genes_in_reg <- regulons[[reg_name]]
  # Filter to only genes in background (tested genes)
  valid_genes <- intersect(genes_in_reg, background)

  if (length(valid_genes) > 0) {
    data.frame(
      TERM = rep(reg_name, length(valid_genes)),
      Gene = valid_genes,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}) %>%
  bind_rows()

cat("[05b_regulon] Created TERM2GENE with", length(unique(regulon_term2gene$TERM)),
    "regulons and", nrow(regulon_term2gene), "gene-regulon associations\n")

# Check for regulon genes not in background
all_regulon_genes <- unique(unlist(regulons))
background_set <- unique(background)
unmapped_regulon_genes <- setdiff(all_regulon_genes, background_set)
if (length(unmapped_regulon_genes) > 0) {
  cat("[05b_regulon] Warning:", length(unmapped_regulon_genes),
      "regulon genes not in tested background\n")
}

# --- Gene sets: same as in 05_enrichment_go_pathway.R ---
gene_sets <- list(
  "RNA-seq up-regulated"    = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  >  logfc_cutoff) %>% pull(name),
  "RNA-seq down-regulated"  = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  < -logfc_cutoff) %>% pull(name),
  "Ribo-seq up-regulated"   = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo >  logfc_cutoff) %>% pull(name),
  "Ribo-seq down-regulated" = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo < -logfc_cutoff) %>% pull(name)
)

# --- Run ORA for each gene set against regulon definitions ---
out_dir <- file.path(dir_enrichment, "regulon")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- lapply(names(gene_sets), function(geneset_name) {
  sig_genes <- gene_sets[[geneset_name]]
  cat("[05b_regulon] ORA for:", geneset_name, "(n =", length(sig_genes), "genes)\n")

  if (length(sig_genes) == 0) {
    cat("[05b_regulon]   Empty gene set, skipping\n")
    return(NULL)
  }

  ora <- enricher(
    gene = sig_genes,
    TERM2GENE = regulon_term2gene,
    universe = background,
    pAdjustMethod = "fdr",
    minGSSize = ora_min_gs_size,
    pvalueCutoff = ora_pvalue_cutoff,
    qvalueCutoff = ora_qvalue_cutoff
  )

  if (is.null(ora) || nrow(as.data.frame(ora)) == 0) {
    cat("[05b_regulon]   No enriched regulons\n")
    return(NULL)
  }

  res_df <- as.data.frame(ora) %>%
    mutate(GeneSet = geneset_name, Annotation = "Regulon")

  # Save TSV
  safe_name <- gsub("[^a-zA-Z0-9_]+", "_", geneset_name)
  write.table(res_df, file.path(out_dir, paste0(safe_name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Save dotplot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- tryCatch(
      clusterProfiler::dotplot(ora, showCategory = 20) +
        ggplot2::ggtitle(paste(geneset_name, "- Regulon")) +
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

regulon_ora_results <- bind_rows(all_results)

# Save combined results
write.table(regulon_ora_results, file.path(dir_enrichment, "regulon_ora_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[05b_regulon] Done. ", nrow(regulon_ora_results), "enriched regulons written to", out_dir, "\n")
