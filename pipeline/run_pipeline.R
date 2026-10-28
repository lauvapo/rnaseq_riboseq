# ============================================================
# Master pipeline driver.
#
# Runs, in order: quality control (replicate correlations & RNA-Ribo
# correlation) -> differential expression -> transcription network topology
# -> GO/pathway ORA (DE-based) -> regulon ORA -> iModulon ORA -> spreadsheet
# export -> summary figures.
#
# Usage:
#   Rscript pipeline/run_pipeline.R
#
# Prerequisite (upstream of this pipeline, not included here): raw
# FASTQ reads already trimmed, aligned, deduplicated, and counted per
# CDS with featureCounts into RNA_seq/counts/CDS and Ribo_seq/counts/CDS.
# ============================================================

# Finds this script's own folder from R's command-line arguments (the
# "--file=..." flag Rscript passes itself), so `source()` below can
# find the other pipeline/*.R scripts by relative path regardless of
# which directory `Rscript pipeline/run_pipeline.R` was launched from.
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
pipeline_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else "pipeline"

# Every step script runs in this same R session (via source()), so
# later steps can see variables/functions defined by earlier ones
# (e.g. 06 reuses `all_info` loaded by 00; 07/09 reuse the ORA helpers
# loaded by 05 via 06). The banner just makes long log output
# skimmable when something fails partway through.
run_step <- function(script_name) {
  cat("\n================================================================\n")
  cat("STEP:", script_name, "\n")
  cat("================================================================\n")
  source(file.path(pipeline_dir, script_name), echo = FALSE)
}

source(file.path(pipeline_dir, "00_config.R"))

run_step("02_quality_control_correlations.R")
run_step("03_differential_expression.R")
run_step("04_upset_plot.R")
run_step("05b_volcano_plots.R")
run_step("05c_scatter_transcriptome_translatome.R")
run_step("05d_gene_category_heatmaps.R")
run_step("05_transcription_network_analysis.R")
run_step("05_enrichment_go_pathway.R")
run_step("05_enrichment_regulon.R")
run_step("05_enrichment_imodulons.R")
run_step("06_export_spreadsheets.R")
run_step("07_summary_figures.R")
run_step("08_enrichment_heatmaps.R")
run_step("09_enrichment_logfc_heatmaps.R")
run_step("10_gene_logfc_heatmaps.R")
run_step("11_curated_gene_heatmap.R")
run_step("12_pathway_ridgeplots.R")

cat("\n================================================================\n")
cat("Pipeline finished!\n")
cat("================================================================\n")
cat("Results workbook:", file.path(dir_results, "DE_and_enrichment_results.xlsx"), "\n")
cat("Summary figures:", file.path(dir_results, "*.png"), "\n")
cat("Enrichment heatmaps:", file.path(dir_results, "*heatmap*.png"), "\n")
