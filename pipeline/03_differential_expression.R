# ============================================================
# Step 3: Differential expression (RNA-seq and Ribo-seq)
#
# Loads CDS count matrices, runs DESeq2 for each assay (no manual
# low-count pre-filter - DESeq2's own independent filtering during
# results() already down-weights/excludes unreliable low-count genes
# as part of its FDR correction), categorizes genes by cross-assay
# behaviour (Co-regulated / Only RNA-seq / Only Ribo-seq / Opposite /
# Not significant), and writes the tables consumed by every downstream
# step:
#   significative_genes/rna_sig.tsv
#   significative_genes/ribo_sig.tsv
#   significative_genes/all_info.tsv   (per-gene RNA+Ribo summary)
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(edgeR)
})

# DESeq2's standard three-stage workflow: DESeqDataSetFromMatrix() packages
# up the raw counts + which sample is which condition; DESeq() then (1)
# estimates size factors to correct for differing sequencing depth
# between samples, (2) estimates each gene's dispersion (how variable
# its counts are beyond Poisson noise, borrowing statistical strength
# across genes with similar expression level), and (3) fits a
# negative-binomial GLM and runs a Wald test for the condition effect.
# results() extracts that test as one row per gene: log2FoldChange,
# p-value, and padj (Benjamini-Hochberg FDR-adjusted p-value).
run_deseq2 <- function(counts_mat, metadata) {
  stopifnot(levels(metadata$condition)[1] == "Control") # reference level, so a positive log2FoldChange means "up in Cap7"
  dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = metadata, design = ~condition)
  dds <- DESeq(dds)
  list(dds = dds, res = DESeq2::results(dds, pAdjustMethod = "fdr"))
}

# Writes the exact raw CDS count matrix handed to DESeq2 for one assay -
# i.e. every gene actually tested, not just the significant subset -
# so "which genes were used and what were their raw counts" is always
# traceable to a concrete file instead of only the DE result tables.
write_gene_counts <- function(counts_mat, basename) {
  out <- data.frame(gene = rownames(counts_mat), gene_symbol = map_to_gene_symbol(rownames(counts_mat)),
                     counts_mat, check.names = FALSE)
  write.table(out, file.path(dir_sig_genes, basename), sep = "\t", quote = FALSE, row.names = FALSE)
}

# --- RNA-seq ---
rna_counts <- read_count_matrix(dir_rna_counts_cds, rna_samples)
rna_metadata_used <- rna_metadata[colnames(rna_counts), , drop = FALSE]
cat("[02_DE] RNA-seq genes tested:", nrow(rna_counts), "\n")
write_gene_counts(rna_counts, "rna_gene_counts.tsv")

rna_res <- run_deseq2(rna_counts, rna_metadata_used)
rna_res_df <- as.data.frame(rna_res$res)
rna_res_df$gene <- rownames(rna_res_df)

rna_sig <- rna_res_df[!is.na(rna_res_df$padj) &
                         abs(rna_res_df$log2FoldChange) > logfc_cutoff &
                         rna_res_df$padj < fdr_cutoff, ]
cat("[02_DE] RNA-seq significant genes:", nrow(rna_sig), "/", nrow(rna_res_df), "\n")
rna_sig$gene_symbol <- map_to_gene_symbol(rownames(rna_sig))
rna_sig <- rna_sig[, c("gene", "gene_symbol", setdiff(colnames(rna_sig), c("gene", "gene_symbol")))]
write.table(rna_sig, file.path(dir_sig_genes, "rna_sig.tsv"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

# --- Ribo-seq ---
ribo_counts <- read_count_matrix(dir_ribo_counts_cds, ribo_samples)
ribo_metadata_used <- ribo_metadata[colnames(ribo_counts), , drop = FALSE]
cat("[02_DE] Ribo-seq genes tested:", nrow(ribo_counts), "\n")
write_gene_counts(ribo_counts, "ribo_gene_counts.tsv")

ribo_res <- run_deseq2(ribo_counts, ribo_metadata_used)
ribo_res_df <- as.data.frame(ribo_res$res)
ribo_res_df$gene <- rownames(ribo_res_df)

ribo_sig <- ribo_res_df[!is.na(ribo_res_df$padj) &
                           abs(ribo_res_df$log2FoldChange) > logfc_cutoff &
                           ribo_res_df$padj < fdr_cutoff, ]
cat("[02_DE] Ribo-seq significant genes:", nrow(ribo_sig), "/", nrow(ribo_res_df), "\n")
ribo_sig$gene_symbol <- map_to_gene_symbol(rownames(ribo_sig))
ribo_sig <- ribo_sig[, c("gene", "gene_symbol", setdiff(colnames(ribo_sig), c("gene", "gene_symbol")))]
write.table(ribo_sig, file.path(dir_sig_genes, "ribo_sig.tsv"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

# --- Spearman correlation between biological replicates (QC) ---
correlation_table <- data.frame(
  Dataset   = c("RNA-seq", "RNA-seq", "Ribo-seq", "Ribo-seq"),
  Condition = c("Control", "Capreomycin", "Control", "Capreomycin"),
  Spearman_correlation = c(
    cor(rna_counts[[rna_control_reps[1]]], rna_counts[[rna_control_reps[2]]], method = "spearman"),
    cor(rna_counts[[rna_cap_reps[1]]],     rna_counts[[rna_cap_reps[2]]],     method = "spearman"),
    cor(ribo_counts[[ribo_control_reps[1]]], ribo_counts[[ribo_control_reps[2]]], method = "spearman"),
    cor(ribo_counts[[ribo_cap_reps[1]]],     ribo_counts[[ribo_cap_reps[2]]],     method = "spearman")
  )
)
correlation_table$Spearman_correlation <- round(correlation_table$Spearman_correlation, 3)
write.table(correlation_table, file.path(dir_results, "replicate_correlation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Cross-assay categorization + all_info.tsv ---
# Bring the two assays' DESeq2 results into one table, one row per gene,
# with RNA-seq and Ribo-seq stats side by side (full_join keeps every
# gene tested in EITHER assay, even if it was filtered out of the
# other - those cells are just NA).
rna_slim  <- rna_res_df[, c("gene", "log2FoldChange", "padj")]
colnames(rna_slim)[2:3] <- c("logFC_rna", "padj_rna")
ribo_slim <- ribo_res_df[, c("gene", "log2FoldChange", "padj")]
colnames(ribo_slim)[2:3] <- c("logFC_ribo", "padj_ribo")

merged_df <- full_join(rna_slim, ribo_slim, by = "gene")
# How differently-sized the two assays' effects are for this gene - near
# 0 means RNA-seq and Ribo-seq logFC are about the same magnitude
# (translation is just tracking transcription); large means one assay
# shows a much bigger swing than the other.
merged_df$extremeness <- abs(abs(merged_df$logFC_rna) - abs(merged_df$logFC_ribo))

sig_rna  <- !is.na(merged_df$padj_rna)  & abs(merged_df$logFC_rna)  > logfc_cutoff & merged_df$padj_rna  < fdr_cutoff
sig_ribo <- !is.na(merged_df$padj_ribo) & abs(merged_df$logFC_ribo) > logfc_cutoff & merged_df$padj_ribo < fdr_cutoff
same_direction <- sign(merged_df$logFC_rna) == sign(merged_df$logFC_ribo)

# Categorize every gene by which assay(s) called it significant and
# whether the direction of change agrees between them. This is a
# descriptive cross-tabulation, not a formal statistical test for
# translational regulation - see 08_translation_efficiency.R for that.
#   Co-regulated:  significant in both assays, same direction (e.g. up in both)
#   Opposite:      significant in both assays, but opposite direction
#   Only RNA-Seq:  significant in RNA-seq only
#   Only Ribo-Seq: significant in Ribo-seq only
merged_df$category <- "Not significant"
merged_df$category[sig_rna & sig_ribo & same_direction]  <- "Co-regulated"
merged_df$category[sig_rna & sig_ribo & !same_direction] <- "Opposite"
merged_df$category[!sig_rna & sig_ribo] <- "Only Ribo-Seq"
merged_df$category[sig_rna & !sig_ribo] <- "Only RNA-Seq"

# Which quadrant of the RNA-vs-Ribo logFC scatter this gene falls in,
# independent of whether either logFC is "significant" - used for
# plotting/labeling downstream, not for the category calls above.
merged_df$direction <- with(merged_df, ifelse(
  logFC_rna > 0 & logFC_ribo > 0, "pos_pos", ifelse(
    logFC_rna < 0 & logFC_ribo < 0, "neg_neg", ifelse(
      logFC_rna > 0 & logFC_ribo < 0, "pos_neg", "neg_pos"))))

all_info <- merged_df
# featureCounts' Geneid (the "gene" column) is a locus tag (e.g.
# "Rv3252c") for every gene in this project's count files - there's no
# actual symbol-to-locus_tag conversion happening here for THIS
# dataset, but map_to_locus_tag is kept for robustness in case a future
# count file ever uses gene symbols as its Geneid instead.
all_info$name <- map_to_locus_tag(all_info$gene, gene_map)
# Human-readable label (e.g. "alkB") for spreadsheets and plots; "name"
# above stays the locus tag and remains the join key used throughout
# the rest of the pipeline (see 00_config.R for why).
all_info$gene_symbol <- map_to_gene_symbol(all_info$name)
all_info <- all_info[, c("gene", "name", "gene_symbol", setdiff(colnames(all_info), c("gene", "name", "gene_symbol")))]
write.table(all_info, file.path(dir_sig_genes, "all_info.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Summary counts table: one-glance gene tallies for the paper/thesis text ---
summary_table <- data.frame(
  Metric = c("Total genes", "Significant genes", "Up", "Down",
             "Significant Up", "Significant Down", "Co-regulated",
             "Opposite (both sig)", "Only RNA-Seq", "Only Ribo-Seq"),
  RNA_seq = c(
    nrow(rna_res_df), nrow(rna_sig),
    sum(rna_res_df$log2FoldChange > 0, na.rm = TRUE),
    sum(rna_res_df$log2FoldChange < 0, na.rm = TRUE),
    sum(rna_sig$log2FoldChange > 0), sum(rna_sig$log2FoldChange < 0),
    sum(merged_df$category == "Co-regulated"), sum(merged_df$category == "Opposite"),
    sum(merged_df$category == "Only RNA-Seq"), NA
  ),
  Ribo_seq = c(
    nrow(ribo_res_df), nrow(ribo_sig),
    sum(ribo_res_df$log2FoldChange > 0, na.rm = TRUE),
    sum(ribo_res_df$log2FoldChange < 0, na.rm = TRUE),
    sum(ribo_sig$log2FoldChange > 0), sum(ribo_sig$log2FoldChange < 0),
    sum(merged_df$category == "Co-regulated"), sum(merged_df$category == "Opposite"),
    NA, sum(merged_df$category == "Only Ribo-Seq")
  )
)
write.table(summary_table, file.path(dir_results, "DE_summary_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[02_DE] Done. all_info.tsv, rna_sig.tsv, ribo_sig.tsv written to", dir_sig_genes, "\n")
