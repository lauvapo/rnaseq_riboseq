# ============================================================
# Pipeline configuration: shared paths, sample design, and
# statistical parameters used by every other pipeline/*.R script.
# Source this file first (run_pipeline.R does this automatically).
#
# Nothing in this file runs an analysis - it only defines WHERE things
# are (file paths), WHAT the experiment looks like (which samples are
# Control vs Cap7), and WHICH thresholds every later step should use,
# so a parameter only ever needs to change in one place.
# ============================================================
# clear the environment to avoid accidental carryover
# Absolute path to the project root. Every other path below is built
# from this one with file.path(), so the whole pipeline is portable to
# a new machine by editing only this single line.
project_dir <- "/home/lauvapo/Documentos/tesis/paper_capreomycin_remapping"

# --- Gene map (copied locally for reproducibility) ---
# Maps gene symbols (e.g. "alkB") to Mtb H37Rv locus tags (e.g. "Rv3252c").
# Originally sourced from an external file shared with other, unrelated
# analyses on this machine (~/Documentos/tesis/riboseq_def/gene_map.tsv);
# this project-local copy additionally has ~143 of its ~145 NA locus_tag
# rows (mostly tRNA symbols and some newer gene names, e.g. the
# ergothioneine/mycofactocin pathway genes) filled in from a fuller
# Mycobrowser v5 export (Mycobacterium_tuberculosis_H37Rv_txt_v5.txt, also
# project-local). Regenerate by re-running that patch against a fresh copy
# of the external file if it's ever updated upstream.
gene_map_file <- file.path(project_dir, "pathways", "gene_map.tsv")

# --- Project-local directories (created below if missing) ---
# Inputs (read-only for this pipeline; produced upstream by trimming +
# alignment + featureCounts, outside these scripts):
dir_rna_counts_cds    <- file.path(project_dir, "data/RNA_seq/counts/CDS")   # per-gene RNA-seq read counts
dir_ribo_counts_cds   <- file.path(project_dir, "data/Ribo_seq/counts/CDS")  # per-gene Ribo-seq read counts
dir_rna_read_lengths  <- file.path(project_dir, "data/RNA_seq/read_lengths")
dir_ribo_read_lengths <- file.path(project_dir, "data/Ribo_seq/read_lengths")
# Outputs (this pipeline writes here; see results/README.md for what's in each file):
dir_sig_genes         <- file.path(project_dir, "results/significative_genes")  # per-gene result tables (DE)
dir_enrichment        <- file.path(project_dir, "results/enrichment_results")   # GO/pathway ORA tables and dotplots
dir_pathways          <- file.path(project_dir, "pathways")             # GO/pathway annotation inputs
dir_results           <- file.path(project_dir, "results")              # the final Excel workbook and summary tables

# All pipeline outputs live under these folders only (never in the
# legacy images/, plots/ folders owned by the standalone .Rmd notebooks).
for (d in c(dir_sig_genes, dir_enrichment,
            file.path(dir_enrichment, "go_pathway"),
            file.path(dir_enrichment, "network"),
            file.path(dir_enrichment, "regulon"),
            file.path(dir_enrichment, "network_enrichment"),
            file.path(dir_enrichment, "imodulons"), dir_results)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- GO term / pathway annotation used for ORA (copied locally for
#     reproducibility; originally from Mycobrowser "All instances of
#     Genes/Pathways in M. tuberculosis H37Rv" exports) ---
go_annotation_file      <- file.path(dir_pathways, "go_terms_genes.txt")
pathway_annotation_file <- file.path(dir_pathways, "pathway_genes.txt")

# --- Regulon and transcription network files for regulon enrichment ---
regulon_annotation_file  <- file.path(dir_pathways, "selected_regulon.gmt")
transcrnetwork_file      <- file.path(dir_pathways, "transcrnetwork_modulome.csv")

# --- Sample design ---
# Two biological replicates per condition, per assay (2 Control + 2
# Cap7 = 4 samples for RNA-seq, and separately 4 for Ribo-seq). Sample
# codes are the library IDs used in every filename in RNA_seq/ and
# Ribo_seq/ (e.g. "J_S4_CDS.dedup.csv").
rna_samples  <- c("J_S4", "M_S4", "L_S6", "O_S6")   # J_S4, M_S4 = Control; L_S6, O_S6 = capreomycin (Cap7)
ribo_samples <- c("D_S1", "G_S1", "F_S3", "I_S3")   # D_S1, G_S1 = Control; F_S3, I_S3 = capreomycin (Cap7)

# "Control" is set as the first factor level here (rather than via a
# post-hoc relevel() call) so DESeqDataSetFromMatrix() picks it up as
# the reference level immediately, with no relevel warning.
condition_levels <- c("Control", "Cap7")

rna_metadata <- data.frame(
  sample = rna_samples,
  condition = factor(ifelse(rna_samples %in% c("J_S4", "M_S4"), "Control", "Cap7"),
                      levels = condition_levels),
  row.names = rna_samples
)

ribo_metadata <- data.frame(
  sample = ribo_samples,
  condition = factor(ifelse(ribo_samples %in% c("D_S1", "G_S1"), "Control", "Cap7"),
                      levels = condition_levels),
  row.names = ribo_samples
)

# Split out again by condition, since almost every filtering/statistics
# step needs "the two Control replicates" and "the two Cap7 replicates"
# as separate groups rather than the full 4-sample vector above.
rna_control_reps  <- c("J_S4", "M_S4")
rna_cap_reps      <- c("L_S6", "O_S6")
ribo_control_reps <- c("D_S1", "G_S1")
ribo_cap_reps     <- c("F_S3", "I_S3")

# --- Statistical thresholds (shared by DE + enrichment steps) ---
# A gene counts as "significant" (DE, translation efficiency, etc.)
# when BOTH conditions hold: adjusted p-value (FDR) below fdr_cutoff,
# AND effect size |log2 fold-change| above logfc_cutoff. Requiring both
# avoids calling a gene "significant" purely because it has a tiny,
# biologically negligible fold-change measured very precisely.
fdr_cutoff   <- 0.1   # 10% false discovery rate
logfc_cutoff <- 0.5   # ~1.4x fold-change in either direction

# --- ORA (over-representation analysis / enrichment) parameters ---
# minGSSize excludes GO/pathway terms with fewer than this many genes
# annotated (too small to test meaningfully). pvalueCutoff/qvalueCutoff
# are the significance thresholds enricher() applies when deciding
# which enriched terms to report.
ora_pvalue_cutoff <- 0.1
ora_qvalue_cutoff <- 0.1
ora_min_gs_size   <- 3

# --- Helper: read one featureCounts CDS count file into a 2-col df ---
# featureCounts output has several metadata columns (Chr, Start, End,
# Strand, Length) before the actual count column, and the count column
# is named after the BAM file it came from (e.g. "./alignments/J_S4...bam")
# rather than the sample name - this keeps only Geneid + the counts,
# and renames the count column to the plain sample name (e.g. "J_S4").
read_one_featurecounts <- function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t", comment.char = "#")
  count_col <- grep("\\.bam$", colnames(df), value = TRUE)
  df2 <- df[, c("Geneid", count_col)]
  colnames(df2)[2] <- tools::file_path_sans_ext(basename(file))
  df2
}

# --- Helper: build a gene x sample count matrix from a set of CDS files ---
# Reads one featureCounts file per sample and joins them side by side
# into a single genes-as-rows, samples-as-columns matrix (the format
# DESeq2 expects as its countData argument). A gene missing from one
# sample's file (shouldn't normally happen, since all samples are
# counted against the same genome annotation) would produce an NA from
# the join - those are set to 0 rather than silently dropped.
read_count_matrix <- function(dir, samples, suffix = "_CDS.dedup.csv") {
  files <- file.path(dir, paste0(samples, suffix))
  stopifnot(all(file.exists(files)))
  counts_list <- lapply(files, read_one_featurecounts)
  counts_all  <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), counts_list)
  rownames(counts_all) <- counts_all$Geneid
  mat <- counts_all[, -1, drop = FALSE]
  mat[] <- lapply(mat, function(x) round(as.numeric(x)))
  mat[is.na(mat)] <- 0
  colnames(mat) <- gsub(paste0(suffix, "$"), "", colnames(mat))
  colnames(mat) <- gsub("_CDS\\.dedup$", "", colnames(mat))
  as.data.frame(mat)
}

# --- Gene map (name <-> locus_tag), loaded once and shared by every step ---
gene_map <- read.delim(gene_map_file, header = TRUE, sep = "\t")

# --- Helper: map gene symbols to locus tags via the project gene map.
#     Genes already given as a locus tag (or with no map entry) pass
#     through unchanged. ---
map_to_locus_tag <- function(genes, gene_map) {
  lookup <- c(setNames(gene_map$locus_tag, gene_map$name),
              setNames(gene_map$locus_tag, gene_map$locus_tag))
  ifelse(genes %in% names(lookup), unname(lookup[genes]), genes)
}

# --- Helper: map locus tags to a human-readable gene SYMBOL, for display
#     only (plots, spreadsheet labels) - never use this for joins.
#
#     gene_map.tsv is many-to-one: each locus tag has several aliases mixed
#     together under "name" (a protein-accession ID like "CCP42739.1", the
#     real gene symbol like "pstP", and/or a self-referential row where
#     name == locus_tag for genes with no assigned symbol at all). We keep
#     only candidates that are neither self-referential nor an accession
#     ID; a handful of locus tags (~1%) still have >1 candidate afterward,
#     resolved by just taking the first. Genes with no real symbol fall
#     back to their locus tag, since every gene has one but not every gene
#     has a symbol. ---
build_symbol_lookup <- function(gene_map) {
  candidates <- gene_map[
    !is.na(gene_map$locus_tag) &
      gene_map$name != gene_map$locus_tag &
      !grepl("^CCP[0-9]+\\.[0-9]+$", gene_map$name),
  ]
  candidates <- candidates[!duplicated(candidates$locus_tag), ]
  setNames(candidates$name, candidates$locus_tag)
}

gene_symbol_lookup <- build_symbol_lookup(gene_map)

map_to_gene_symbol <- function(locus_tags) {
  symbol <- unname(gene_symbol_lookup[locus_tags])
  ifelse(is.na(symbol), locus_tags, symbol)
}

cat("[00_config] Configuration loaded. Project dir:", project_dir, "\n")
