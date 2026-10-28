# ============================================================
# Step 5: GO / pathway ORA for the DESeq2-significant gene sets
#
# Over-representation analysis (GO Biological Process / Cellular
# Component / Molecular Function + custom pathways) for the classic
# DESeq2 up/down-regulated gene sets (RNA-seq and Ribo-seq).
#
# Outputs one TSV (+ dotplot PNG) per gene-set x annotation
# combination under enrichment_results/, plus a combined table:
#   enrichment_results/go_pathway_ora_results.tsv
# ============================================================

source(file.path(pipeline_dir, "04_enrichment_common.R"))
# (defines `all_info`, `background`, and the annotation/ORA helpers used below -
# see that file for what ORA/hypergeometric enrichment actually computes)

# Gene sets = genes significant in that assay's own DESeq2 test, split
# by direction (independent of whether the other assay also called
# them significant, so Co-regulated genes are included here too, same
# as the original analysis).
gene_sets <- list(
  "RNA-seq up-regulated"    = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  >  logfc_cutoff) %>% pull(name),
  "RNA-seq down-regulated"  = all_info %>% filter(!is.na(padj_rna),  padj_rna  < fdr_cutoff, logFC_rna  < -logfc_cutoff) %>% pull(name),
  "Ribo-seq up-regulated"   = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo >  logfc_cutoff) %>% pull(name),
  "Ribo-seq down-regulated" = all_info %>% filter(!is.na(padj_ribo), padj_ribo < fdr_cutoff, logFC_ribo < -logfc_cutoff) %>% pull(name)
)

out_dir <- file.path(dir_enrichment, "go_pathway")
all_results <- lapply(names(gene_sets), function(name) {
  cat("[05_enrichment] ORA for:", name, "(n =", length(gene_sets[[name]]), "genes)\n")
  run_ora_all_annotations(gene_sets[[name]], name, out_dir)
})

go_pathway_ora_results <- bind_rows(all_results)
write.table(go_pathway_ora_results, file.path(dir_enrichment, "go_pathway_ora_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[05_enrichment] Done. ", nrow(go_pathway_ora_results), "enriched terms written to", out_dir, "\n")
