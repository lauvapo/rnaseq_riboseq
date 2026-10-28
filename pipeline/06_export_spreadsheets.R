# ============================================================
# Step 11: Export DE + enrichment results into one Excel workbook
#
# Gathers every result table produced by steps 01-10 into a single
# multi-sheet workbook:
#   results/DE_and_enrichment_results.xlsx
# ============================================================

suppressPackageStartupMessages(library(openxlsx))

# Any given result file might legitimately not exist (e.g. an ORA step
# found nothing significant anywhere, so it never wrote an output file)
# - this reads a file if present, or returns NULL to skip it cleanly
# instead of erroring, since a missing file here is an expected outcome,
# not a bug.
read_tsv_safe <- function(path, ...) {
  if (!file.exists(path)) {
    message("[11_export] Missing (skipped): ", path)
    return(NULL)
  }
  read.delim(path, ...)
}

# --- ORA result tables get split by:
#     1. Annotation type (GO_BP, GO_CC, GO_MF, Pathway, Regulon, Network, iModulon)
#     2. GeneSet column contains "RNA-seq up/down" or "Ribo-seq up/down"
annotation_sheet_names <- c(
  "Biological Process" = "GO_BP",
  "Cellular Component" = "GO_CC",
  "Molecular Function" = "GO_MF",
  "Pathway"             = "Pathway",
  "Regulon"            = "Regulon",
  "iModulon"           = "iModulon"
)

# Splits ORA results by annotation type and GeneSet (which contains assay + direction)
split_ora_by_annotation_and_geneset <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(list())
  df <- df[order(df$GeneSet, df$p.adjust), ]

  result_list <- list()

  # Split by annotation type
  by_annotation <- split(df, df$Annotation, drop = TRUE)

  for (annot in names(by_annotation)) {
    annot_data <- by_annotation[[annot]]
    annot_code <- annotation_sheet_names[annot]

    # Split by GeneSet (e.g., "RNA-seq up-regulated", "Ribo-seq down-regulated")
    by_geneset <- split(annot_data, annot_data$GeneSet, drop = TRUE)

    for (geneset in names(by_geneset)) {
      geneset_data <- by_geneset[[geneset]]
      # Convert "RNA-seq up-regulated" -> "RNA_up", "Ribo-seq down-regulated" -> "Ribo_down"
      geneset_short <- gsub("seq ", "", geneset)  # "RNA-seq" -> "RNA", "Ribo-seq" -> "Ribo"
      geneset_short <- gsub("-regulated", "", geneset_short)  # "up-regulated" -> "up"
      geneset_short <- tolower(gsub("-", "_", geneset_short))  # "RNA up" -> "rna_up"

      sheet_name <- paste0(annot_code, "_", geneset_short)
      result_list[[sheet_name]] <- geneset_data
    }
  }

  result_list
}

ora_sheets <- split_ora_by_annotation_and_geneset(
  read_tsv_safe(file.path(dir_enrichment, "go_pathway_ora_results.tsv"))
)

# --- Regulon enrichment sheets ---
regulon_sheets <- split_ora_by_annotation_and_geneset(
  read_tsv_safe(file.path(dir_enrichment, "regulon_ora_results.tsv"))
)


# --- iModulon enrichment sheets ---
imodulons_sheets <- split_ora_by_annotation_and_geneset(
  read_tsv_safe(file.path(dir_enrichment, "imodulons_ora_results.tsv"))
)

sheets <- c(
  list(
    "All_DE"              = read_tsv_safe(file.path(dir_sig_genes, "all_info.tsv")),
    "DE_summary"          = read_tsv_safe(file.path(dir_results, "DE_summary_counts.tsv")),
    "Replicate_correlation" = read_tsv_safe(file.path(dir_results, "replicate_correlation.tsv")),
    "Read_length_summary" = read_tsv_safe(file.path(dir_results, "read_length_summary.tsv"))
  ),
  ora_sheets,
  regulon_sheets,
  imodulons_sheets
)

sheets <- Filter(function(x) !is.null(x) && nrow(x) > 0, sheets)

wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "Bottom")

# One sheet per entry in `sheets`, with a bold header row, the header
# row frozen so it stays visible while scrolling, auto-sized columns,
# and an Excel autofilter dropdown on every column - so each sheet is
# independently sortable/filterable without needing separate sheets for
# every possible slice of the data.
for (sheet_name in names(sheets)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, sheets[[sheet_name]], headerStyle = header_style)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(sheets[[sheet_name]])), widths = "auto")
  addFilter(wb, sheet_name, rows = 1, cols = seq_len(ncol(sheets[[sheet_name]])))
}

out_file <- file.path(dir_results, "DE_and_enrichment_results.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)

cat("[11_export] Workbook written:", out_file, "\n")
cat("[11_export] Sheets included:", paste(names(sheets), collapse = ", "), "\n")
