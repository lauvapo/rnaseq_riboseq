# ============================================================
# Step 2: Quality Control - Replicate & RNA-Ribo Correlations
#
# Assesses data quality by:
#   1. Computing Spearman correlations between biological
#      replicates within each condition for both RNA-seq and
#      Ribo-seq data
#   2. Comparing transcript abundance (RNA-seq) vs ribosome
#      occupancy (Ribo-seq) across genes by condition
#
# Output files:
#   - qc_replicate_correlations.csv (summary table)
#   - qc_correlation_plot.png (replicate correlations)
#   - qc_correlation_plot.svg (replicate correlations, vector)
#   - rna_ribo_correlation_scatter.png (RNA vs Ribo scatter)
#   - rna_ribo_correlation_scatter.svg (RNA vs Ribo, vector)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# Helper function to load counts from a featureCounts output file
# (already defined in 00_config.R; the guard just avoids clobbering it)
if (!exists("read_one_featurecounts")) {
  read_one_featurecounts <- function(file) {
    lines <- readLines(file)
    data <- read.delim(text = lines, stringsAsFactors = FALSE, row.names = 1)
    list(rownames(data), data[[1]])
  }
}

# Load RNA-seq counts
cat("[04_qc_correlations] Loading RNA-seq data...\n")
rna_counts <- list()
for (sample in rna_samples) {
  file <- file.path(dir_rna_counts_cds, paste0(sample, "_CDS.dedup.csv"))
  if (!file.exists(file)) {
    warning("RNA file not found: ", file)
    next
  }
  data <- read_one_featurecounts(file)
  rna_counts[[sample]] <- setNames(data[[2]], data[[1]])
}

# Load Ribo-seq counts
cat("[04_qc_correlations] Loading Ribo-seq data...\n")
ribo_counts <- list()
for (sample in ribo_samples) {
  file <- file.path(dir_ribo_counts_cds, paste0(sample, "_CDS.dedup.csv"))
  if (!file.exists(file)) {
    warning("Ribo file not found: ", file)
    next
  }
  data <- read_one_featurecounts(file)
  ribo_counts[[sample]] <- setNames(data[[2]], data[[1]])
}

# ============================================================
# Calculate correlations for RNA-seq
# ============================================================

rna_spearman <- data.frame(
  comparison = c("Control", "Cap7"),
  correlation = c(
    cor(rna_counts[[rna_control_reps[1]]],
        rna_counts[[rna_control_reps[2]]],
        method = "spearman", use = "complete.obs"),

    cor(rna_counts[[rna_cap_reps[1]]],
        rna_counts[[rna_cap_reps[2]]],
        method = "spearman", use = "complete.obs")
  ),
  dataset = "RNA-seq"
)

# ============================================================
# Calculate correlations for Ribo-seq
# ============================================================

ribo_spearman <- data.frame(
  comparison = c("Control", "Cap7"),
  correlation = c(
    cor(ribo_counts[[ribo_control_reps[1]]],
        ribo_counts[[ribo_control_reps[2]]],
        method = "spearman", use = "complete.obs"),

    cor(ribo_counts[[ribo_cap_reps[1]]],
        ribo_counts[[ribo_cap_reps[2]]],
        method = "spearman", use = "complete.obs")
  ),
  dataset = "Ribo-seq"
)

# ============================================================
# Combine into single dataframe and export
# ============================================================

all_spearman <- rbind(rna_spearman, ribo_spearman)

cat("[04_qc_correlations] Spearman correlations:\n")
print(all_spearman)

# Create a more detailed summary table with sample names
rna_table <- data.frame(
  Dataset = "RNA-seq",
  Condition = c("Control", "Cap7"),
  Replicate_1 = c(rna_control_reps[1], rna_cap_reps[1]),
  Replicate_2 = c(rna_control_reps[2], rna_cap_reps[2]),
  Spearman_correlation = rna_spearman$correlation
)

ribo_table <- data.frame(
  Dataset = "Ribo-seq",
  Condition = c("Control", "Cap7"),
  Replicate_1 = c(ribo_control_reps[1], ribo_cap_reps[1]),
  Replicate_2 = c(ribo_control_reps[2], ribo_cap_reps[2]),
  Spearman_correlation = ribo_spearman$correlation
)

correlation_table <- rbind(rna_table, ribo_table)
correlation_table$Spearman_correlation <- round(correlation_table$Spearman_correlation, 4)

# Export to CSV
write.csv(correlation_table,
          file = file.path(dir_results, "qc_replicate_correlations.csv"),
          row.names = FALSE, quote = FALSE)
cat("[04_qc_correlations] Exported correlation table to: qc_replicate_correlations.csv\n")

# ============================================================
# Generate visualization
# ============================================================

if (!require("ggplot2", quietly = TRUE)) {
  cat("Warning: ggplot2 not found, skipping correlation plot\n")
} else {
  all_spearman$comparison <- factor(all_spearman$comparison,
                                     levels = c("Control", "Cap7"))

  p <- ggplot(all_spearman,
              aes(x = dataset, y = correlation, color = comparison)) +
    geom_point(size = 6, position = position_dodge(width = 0.4)) +
    geom_text(
      aes(label = round(correlation, 4)),
      position = position_dodge(width = 0.4),
      vjust = -1,
      size = 6,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      "Control" = "#0D0887FF",
      "Cap7" = "#E16462FF"
    )) +
    ylim(0, 1) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) +
    labs(
      title = "Spearman correlation between biological replicates",
      x = "",
      y = "Spearman correlation coefficient",
      color = "Condition"
    )

  # Export as PNG
  png(file.path(dir_results, "qc_correlation_plot.png"),
      width = 3000, height = 2500, res = 300)
  print(p)
  dev.off()
  cat("[04_qc_correlations] Exported PNG plot to: qc_correlation_plot.png\n")

  # Export as SVG
  svg(file.path(dir_results, "qc_correlation_plot.svg"),
      width = 10, height = 8)
  print(p)
  dev.off()
  cat("[04_qc_correlations] Exported SVG plot to: qc_correlation_plot.svg\n")
}

cat("[04_qc_correlations] Quality control step complete.\n")

# ============================================================
# RNA-Ribo Correlation Analysis
# ============================================================

cat("[04_qc_correlations] Starting RNA-Ribo correlation analysis...\n")

# Load edgeR for CPM normalization
if (!require("edgeR", quietly = TRUE)) {
  cat("Warning: edgeR not installed. Attempting alternative CPM calculation.\n")
  cpm <- function(x) sweep(x, 2, colSums(x) / 1e6, "/")
}

# Reload counts as matrices for CPM normalization
cat("[04_qc_correlations] Preparing counts for RNA-Ribo correlation...\n")
rna_counts_matrix <- do.call(cbind, rna_counts)
ribo_counts_matrix <- do.call(cbind, ribo_counts)

# Normalize to CPM and log2 transform
cat("[04_qc_correlations] Normalizing counts to log2 CPM...\n")
rna_cpm <- cpm(as.matrix(rna_counts_matrix))
ribo_cpm <- cpm(as.matrix(ribo_counts_matrix))
rna_cpm <- log2(rna_cpm + 1)
ribo_cpm <- log2(ribo_cpm + 1)

# Keep only genes present in both assays
common_genes <- intersect(rownames(rna_cpm), rownames(ribo_cpm))
cat("[04_qc_correlations] Using", length(common_genes),
    "genes present in both RNA-seq and Ribo-seq\n")

rna_cpm <- rna_cpm[common_genes, ]
ribo_cpm <- ribo_cpm[common_genes, ]

# Calculate per-condition averages
rna_ctrl <- rowMeans(rna_cpm[, rna_control_reps])
ribo_ctrl <- rowMeans(ribo_cpm[, ribo_control_reps])

rna_cap <- rowMeans(rna_cpm[, rna_cap_reps])
ribo_cap <- rowMeans(ribo_cpm[, ribo_cap_reps])

# Combine into dataframes for plotting
df_ctrl <- data.frame(
  RNA = rna_ctrl,
  Ribo = ribo_ctrl,
  Condition = "Control"
)

df_cap <- data.frame(
  RNA = rna_cap,
  Ribo = ribo_cap,
  Condition = "Cap7"
)

df_corr <- rbind(df_ctrl, df_cap)
df_corr$Condition <- factor(df_corr$Condition, levels = c("Control", "Cap7"))

# Generate visualization
if (!require("ggplot2", quietly = TRUE)) {
  cat("Warning: ggplot2 not found, skipping RNA-Ribo correlation plot\n")
} else {
  if (!require("ggpubr", quietly = TRUE)) {
    cat("Warning: ggpubr not found. Spearman correlations will not be displayed.\n")
    corr_scatter <- ggplot(df_corr, aes(x = RNA, y = Ribo, color = Condition)) +
      geom_point(alpha = 0.15, size = 0.9) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 1.5) +
      scale_color_manual(values = c(
        "Control" = "#0D0887FF",
        "Cap7" = "#E16462FF"
      )) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")
      ) +
      labs(
        title = "RNA-seq vs Ribo-seq correlation by condition",
        x = "RNA-seq log2 CPM",
        y = "Ribo-seq log2 CPM",
        color = "Condition"
      )
  } else {
    corr_scatter <- ggplot(df_corr, aes(x = RNA, y = Ribo, color = Condition)) +
      geom_point(alpha = 0.15, size = 0.9) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 1.5) +
      ggpubr::stat_cor(
        aes(color = Condition),
        method = "spearman",
        show.legend = FALSE,
        size = 5,
        fontface = "bold"
      ) +
      scale_color_manual(values = c(
        "Control" = "#0D0887FF",
        "Cap7" = "#E16462FF"
      )) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")
      ) +
      labs(
        title = "RNA-seq vs Ribo-seq correlation by condition",
        x = "RNA-seq log2 CPM",
        y = "Ribo-seq log2 CPM",
        color = "Condition"
      )
  }

  # Export as PNG
  ggsave(
    file.path(dir_results, "rna_ribo_correlation_scatter.png"),
    plot = corr_scatter,
    width = 10, height = 6, dpi = 300
  )
  cat("[04_qc_correlations] Exported PNG: rna_ribo_correlation_scatter.png\n")

  # Export as SVG
  ggsave(
    file.path(dir_results, "rna_ribo_correlation_scatter.svg"),
    plot = corr_scatter,
    width = 10, height = 6
  )
  cat("[04_qc_correlations] Exported SVG: rna_ribo_correlation_scatter.svg\n")
}

cat("[04_qc_correlations] Quality control analysis complete.\n")

# ============================================================
# Read Mapping by Category (RNA-seq and Ribo-seq)
# ============================================================

cat("[04_qc_correlations] Analyzing reads by mapping category...\n")

dir_rna_counts_def <- file.path(project_dir, "data/RNA_seq/counts/counts_def")
dir_ribo_counts_def <- file.path(project_dir, "data/Ribo_seq/counts/counts_def")

# RNA-seq categories (no IGR for RNA)
rna_categories <- c("CDS", "5UTR", "ncRNA")
# Ribo-seq categories (includes IGR)
ribo_categories <- c("CDS", "5UTR", "IGR", "ncRNA")

# Helper function to read category read counts
read_category_counts <- function(samples, categories, counts_dir, assay_name) {
  read_data <- data.frame()

  for (sample in samples) {
    for (category in categories) {
      summary_file <- file.path(counts_dir,
                                paste0(sample, "_", category, ".dedup.csv.summary"))
      if (file.exists(summary_file)) {
        # Read summary file (first column is Status, second is count)
        summary_data <- read.delim(summary_file, stringsAsFactors = FALSE, row.names = 1)
        colnames(summary_data) <- "count"

        # Extract Assigned reads
        if ("Assigned" %in% rownames(summary_data)) {
          assigned_reads <- summary_data["Assigned", "count"]

          read_data <- rbind(read_data, data.frame(
            sample = sample,
            category = category,
            assigned_reads = assigned_reads,
            assay = assay_name,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  read_data
}

# ============================================================
# RNA-seq Reads by Category
# ============================================================

cat("[04_qc_correlations] Analyzing RNA-seq reads by mapping category...\n")
rna_read_counts <- read_category_counts(rna_samples, rna_categories,
                                         dir_rna_counts_def, "RNA-seq")

# Map samples to conditions
rna_read_counts$condition <- ifelse(rna_read_counts$sample %in% rna_control_reps,
                                    "Control", "Cap7")
rna_read_counts$condition <- factor(rna_read_counts$condition, levels = c("Control", "Cap7"))
rna_read_counts$category <- factor(rna_read_counts$category,
                                   levels = c("CDS", "5UTR", "ncRNA"))

cat("[04_qc_correlations] RNA-seq reads by category:\n")
print(rna_read_counts)

# ============================================================
# Ribo-seq Reads by Category
# ============================================================

cat("[04_qc_correlations] Analyzing Ribo-seq reads by mapping category...\n")
ribo_read_counts <- read_category_counts(ribo_samples, ribo_categories,
                                          dir_ribo_counts_def, "Ribo-seq")

# Map samples to conditions
ribo_read_counts$condition <- ifelse(ribo_read_counts$sample %in% ribo_control_reps,
                                     "Control", "Cap7")
ribo_read_counts$condition <- factor(ribo_read_counts$condition, levels = c("Control", "Cap7"))
ribo_read_counts$category <- factor(ribo_read_counts$category,
                                    levels = c("CDS", "5UTR", "IGR", "ncRNA"))

cat("[04_qc_correlations] Ribo-seq reads by category:\n")
print(ribo_read_counts)

# ============================================================
# Generate Boxplots
# ============================================================

if (!require("ggplot2", quietly = TRUE)) {
  cat("Warning: ggplot2 not found, skipping category boxplots\n")
} else {
  # RNA-seq boxplot
  p_rna_cat <- ggplot(rna_read_counts, aes(x = category, y = assigned_reads,
                                           fill = condition)) +
    geom_boxplot(alpha = 0.8, outlier.size = 2) +
    geom_point(aes(color = condition), position = position_dodge(width = 0.75),
               size = 3, alpha = 0.6) +
    scale_fill_manual(values = c(
      "Control" = "#0D0887FF",
      "Cap7" = "#E16462FF"
    )) +
    scale_color_manual(values = c(
      "Control" = "#0D0887FF",
      "Cap7" = "#E16462FF"
    )) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    ) +
    labs(
      title = "RNA-seq reads mapping by feature category",
      x = "Feature Category",
      y = "Number of assigned reads",
      fill = "Condition",
      color = "Condition"
    ) +
    scale_y_continuous(labels = scales::comma)

  # Export RNA-seq plots
  ggsave(
    file.path(dir_results, "rna_reads_by_category_boxplot.png"),
    plot = p_rna_cat,
    width = 10, height = 7, dpi = 300
  )
  cat("[04_qc_correlations] Exported PNG: rna_reads_by_category_boxplot.png\n")

  ggsave(
    file.path(dir_results, "rna_reads_by_category_boxplot.svg"),
    plot = p_rna_cat,
    width = 10, height = 7
  )
  cat("[04_qc_correlations] Exported SVG: rna_reads_by_category_boxplot.svg\n")

  # Ribo-seq boxplot
  p_ribo_cat <- ggplot(ribo_read_counts, aes(x = category, y = assigned_reads,
                                             fill = condition)) +
    geom_boxplot(alpha = 0.8, outlier.size = 2) +
    geom_point(aes(color = condition), position = position_dodge(width = 0.75),
               size = 3, alpha = 0.6) +
    scale_fill_manual(values = c(
      "Control" = "#0D0887FF",
      "Cap7" = "#E16462FF"
    )) +
    scale_color_manual(values = c(
      "Control" = "#0D0887FF",
      "Cap7" = "#E16462FF"
    )) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    ) +
    labs(
      title = "Ribo-seq reads mapping by feature category",
      x = "Feature Category",
      y = "Number of assigned reads",
      fill = "Condition",
      color = "Condition"
    ) +
    scale_y_continuous(labels = scales::comma)

  # Export Ribo-seq plots
  ggsave(
    file.path(dir_results, "ribo_reads_by_category_boxplot.png"),
    plot = p_ribo_cat,
    width = 10, height = 7, dpi = 300
  )
  cat("[04_qc_correlations] Exported PNG: ribo_reads_by_category_boxplot.png\n")

  ggsave(
    file.path(dir_results, "ribo_reads_by_category_boxplot.svg"),
    plot = p_ribo_cat,
    width = 10, height = 7
  )
  cat("[04_qc_correlations] Exported SVG: ribo_reads_by_category_boxplot.svg\n")

  # ============================================================
  # Stacked Bar Chart: Fraction of reads by category (comparable)
  # ============================================================

  cat("[04_qc_correlations] Creating comparable fraction plot...\n")

  # Calculate total reads per sample per assay
  rna_totals <- aggregate(assigned_reads ~ sample, rna_read_counts, sum)
  ribo_totals <- aggregate(assigned_reads ~ sample, ribo_read_counts, sum)

  # Calculate fractions
  rna_read_counts$fraction <- NA
  ribo_read_counts$fraction <- NA

  for (i in seq_len(nrow(rna_read_counts))) {
    sample <- rna_read_counts$sample[i]
    total <- rna_totals$assigned_reads[rna_totals$sample == sample]
    rna_read_counts$fraction[i] <- (rna_read_counts$assigned_reads[i] / total) * 100
  }

  for (i in seq_len(nrow(ribo_read_counts))) {
    sample <- ribo_read_counts$sample[i]
    total <- ribo_totals$assigned_reads[ribo_totals$sample == sample]
    ribo_read_counts$fraction[i] <- (ribo_read_counts$assigned_reads[i] / total) * 100
  }

  # Combine and average by assay and category
  all_fractions <- rbind(rna_read_counts, ribo_read_counts)

  frac_summary <- all_fractions %>%
    group_by(assay, category) %>%
    summarise(
      mean_fraction = mean(fraction, na.rm = TRUE),
      sd_fraction = sd(fraction, na.rm = TRUE),
      .groups = "drop"
    )

  # Reorder categories for consistent visualization (excluding ncRNA for RNA-seq)
  frac_summary$category <- factor(frac_summary$category,
                                   levels = c("CDS", "5UTR", "IGR", "ncRNA"))
  frac_summary <- frac_summary[order(frac_summary$assay, frac_summary$category), ]

  cat("[04_qc_correlations] Read fractions by feature:\n")
  print(frac_summary)

  # Create grouped (dodged) bar chart
  frac_summary$assay <- factor(frac_summary$assay, levels = c("RNA-seq", "Ribo-seq"))
  frac_summary$category <- factor(frac_summary$category, levels = c("CDS", "5UTR", "ncRNA", "IGR"))

  p_grouped <- ggplot(frac_summary, aes(x = category, y = mean_fraction, fill = assay)) +
    geom_bar(stat = "identity", alpha = 0.85, color = "black", linewidth = 0.8,
             position = position_dodge(width = 0.7)) +
    geom_errorbar(aes(ymin = mean_fraction - sd_fraction,
                      ymax = mean_fraction + sd_fraction),
                  width = 0.2, linewidth = 0.8,
                  position = position_dodge(width = 0.7)) +
    geom_text(aes(label = paste0(round(mean_fraction, 1), "%")),
              position = position_dodge(width = 0.7),
              vjust = -0.5, fontface = "bold", size = 6, color = "black") +
    scale_fill_manual(values = c(
      "RNA-seq" = "#E16462FF",
      "Ribo-seq" = "#0D0887FF"
    )) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
    ) +
    labs(
      title = "Fraction of reads mapping by feature category",
      y = "Fraction of reads (%)",
      fill = "Assay",
      caption = "Note: ncRNA = true small regulatory RNAs only (tRNA/rRNA excluded from annotation)"
    ) +
    ylim(0, 100)

  # Export grouped bar plot
  ggsave(
    file.path(dir_results, "read_fractions_grouped_comparison.png"),
    plot = p_grouped,
    width = 12, height = 7, dpi = 300
  )
  cat("[04_qc_correlations] Exported PNG: read_fractions_grouped_comparison.png\n")

  ggsave(
    file.path(dir_results, "read_fractions_grouped_comparison.svg"),
    plot = p_grouped,
    width = 12, height = 7
  )
  cat("[04_qc_correlations] Exported SVG: read_fractions_grouped_comparison.svg\n")
}

# ============================================================
# Read Length Distribution Analysis
# ============================================================

cat("[04_qc_correlations] Analyzing read length distributions...\n")

dir_rna_read_lengths_def <- file.path(project_dir, "data/RNA_seq/read_lengths/read_lenghts_def")
dir_ribo_read_lengths_def <- file.path(project_dir, "data/Ribo_seq/read_lengths/reads_length_def")

# Helper function to load read length distribution
load_read_lengths <- function(sample, dir) {
  file <- file.path(dir, paste0(sample, ".sortname_mappedread-length-distribution.txt"))
  if (!file.exists(file)) {
    warning("Read length file not found: ", file)
    return(NULL)
  }
  # File format: space-separated "count length", read with read.table
  data <- read.table(file, header = FALSE, col.names = c("count", "length"),
                     stringsAsFactors = FALSE)
  data$count <- as.numeric(data$count)
  data$length <- as.numeric(data$length)
  data$sample <- sample
  data
}

# Load RNA-seq read lengths
cat("[04_qc_correlations] Loading RNA-seq read length distributions...\n")
rna_lengths_list <- lapply(rna_samples, load_read_lengths, dir = dir_rna_read_lengths_def)
rna_lengths_df <- do.call(rbind, Filter(Negate(is.null), rna_lengths_list))
rna_lengths_df$assay <- "RNA-seq"

# Map samples to conditions
rna_lengths_df$condition <- ifelse(rna_lengths_df$sample %in% rna_control_reps,
                                    "Control", "Cap7")
rna_lengths_df$condition <- factor(rna_lengths_df$condition, levels = c("Control", "Cap7"))

# Load Ribo-seq read lengths
cat("[04_qc_correlations] Loading Ribo-seq read length distributions...\n")
ribo_lengths_list <- lapply(ribo_samples, load_read_lengths, dir = dir_ribo_read_lengths_def)
ribo_lengths_df <- do.call(rbind, Filter(Negate(is.null), ribo_lengths_list))
ribo_lengths_df$assay <- "Ribo-seq"

# Map samples to conditions
ribo_lengths_df$condition <- ifelse(ribo_lengths_df$sample %in% ribo_control_reps,
                                     "Control", "Cap7")
ribo_lengths_df$condition <- factor(ribo_lengths_df$condition, levels = c("Control", "Cap7"))

# Combine both
all_lengths_df <- rbind(rna_lengths_df, ribo_lengths_df)

# Calculate summary statistics (without replicating massive arrays)
length_summary <- all_lengths_df %>%
  group_by(assay, sample, condition) %>%
  summarise(
    mean_length = weighted.mean(length, count, na.rm = TRUE),
    total_reads = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

# Add SD for error bars
for (i in seq_len(nrow(length_summary))) {
  subset_data <- all_lengths_df[all_lengths_df$sample == length_summary$sample[i], ]
  mean_len <- length_summary$mean_length[i]
  sd_len <- sqrt(sum(subset_data$count * (subset_data$length - mean_len)^2, na.rm = TRUE) /
                 sum(subset_data$count, na.rm = TRUE))
  length_summary$sd_length[i] <- sd_len
}

cat("[04_qc_correlations] Read length summary statistics:\n")
print(length_summary)

# Save read length summary for Excel export
write.table(length_summary, file = file.path(dir_results, "read_length_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("[04_qc_correlations] Exported read length summary to: read_length_summary.tsv\n")

if (!require("ggplot2", quietly = TRUE)) {
  cat("Warning: ggplot2 not found, skipping read length plots\n")
} else {
  # ============================================================
  # Plot 1: Density plot of read lengths (by assay and condition)
  # ============================================================

  p_density <- ggplot(all_lengths_df, aes(x = length, weight = count, fill = condition)) +
    geom_density(alpha = 0.6, adjust = 1) +
    facet_wrap(~assay, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(
      "Control" = "#0D0887FF",
      "Cap7" = "#E16462FF"
    )) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = "Read length distribution by assay and condition",
      x = "Read length (nt)",
      y = "Density",
      fill = "Condition"
    )

  ggsave(
    file.path(dir_results, "read_length_density_distribution.png"),
    plot = p_density,
    width = 14, height = 6, dpi = 300
  )
  cat("[04_qc_correlations] Exported PNG: read_length_density_distribution.png\n")

  ggsave(
    file.path(dir_results, "read_length_density_distribution.svg"),
    plot = p_density,
    width = 14, height = 6
  )
  cat("[04_qc_correlations] Exported SVG: read_length_density_distribution.svg\n")

  # ============================================================
  # Plot 2: Density plot overlaid by condition
  # ============================================================

  p_violin <- ggplot(all_lengths_df, aes(x = length, weight = count, color = condition,
                                         fill = condition)) +
    geom_density(alpha = 0.3, linewidth = 1.2) +
    facet_wrap(~assay, scales = "free_y", nrow = 1) +
    scale_color_manual(values = c(
      "Control" = "#0D0887FF",
      "Cap7" = "#E16462FF"
    )) +
    scale_fill_manual(values = c(
      "Control" = "#0D0887FF",
      "Cap7" = "#E16462FF"
    )) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = "Read length distribution by condition",
      x = "Read length (nt)",
      y = "Density",
      color = "Condition",
      fill = "Condition"
    )

  ggsave(
    file.path(dir_results, "read_length_condition_comparison.png"),
    plot = p_violin,
    width = 14, height = 6, dpi = 300
  )
  cat("[04_qc_correlations] Exported PNG: read_length_condition_comparison.png\n")

  ggsave(
    file.path(dir_results, "read_length_condition_comparison.svg"),
    plot = p_violin,
    width = 14, height = 6
  )
  cat("[04_qc_correlations] Exported SVG: read_length_condition_comparison.svg\n")

  # ============================================================
  # Plot 3: Mean read length by sample
  # ============================================================

  length_summary$assay_condition <- paste0(length_summary$assay, " - ", length_summary$condition)

  p_barplot <- ggplot(length_summary, aes(x = reorder(sample, mean_length), y = mean_length,
                                          fill = condition)) +
    geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean_length - sd_length, ymax = mean_length + sd_length),
                  width = 0.3, linewidth = 1) +
    facet_wrap(~assay, scales = "free_x") +
    scale_fill_manual(values = c(
      "Control" = "#0D0887FF",
      "Cap7" = "#E16462FF"
    )) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "Mean read length by sample",
      x = "Sample",
      y = "Mean read length (nt)",
      fill = "Condition"
    )

  ggsave(
    file.path(dir_results, "read_length_by_sample.png"),
    plot = p_barplot,
    width = 12, height = 6, dpi = 300
  )
  cat("[04_qc_correlations] Exported PNG: read_length_by_sample.png\n")

  ggsave(
    file.path(dir_results, "read_length_by_sample.svg"),
    plot = p_barplot,
    width = 12, height = 6
  )
  cat("[04_qc_correlations] Exported SVG: read_length_by_sample.svg\n")
}

cat("[04_qc_correlations] Read length analysis complete.\n")

# ============================================================
# Ribo-seq CDS Reads Reproducibility with Read Length Context
# ============================================================

cat("[04_qc_correlations] Analyzing Ribo-seq CDS read reproducibility...\n")

# Load Ribo-seq CDS counts and read length data for correlation analysis
ribo_cds_counts <- list()
for (sample in ribo_samples) {
  file <- file.path(dir_ribo_counts_cds, paste0(sample, "_CDS.dedup.csv"))
  if (file.exists(file)) {
    data <- read_one_featurecounts(file)
    ribo_cds_counts[[sample]] <- setNames(data[[2]], data[[1]])
  }
}

# Get mean read lengths per gene per sample
ribo_read_length_means <- list()
for (sample in ribo_samples) {
  file <- file.path(dir_ribo_read_lengths_def, paste0(sample, ".sortname_mappedread-length-distribution.txt"))
  if (file.exists(file)) {
    lengths <- read.table(file, col.names = c("count", "length"))
    # Calculate weighted mean read length
    weighted_mean_length <- weighted.mean(lengths$length, lengths$count)
    ribo_read_length_means[[sample]] <- weighted_mean_length
  }
}

# Create comparison dataframes for each condition
# Control: D_S1 vs G_S1
if (!is.null(ribo_cds_counts[["D_S1"]]) && !is.null(ribo_cds_counts[["G_S1"]])) {
  control_comparison <- data.frame(
    gene = names(ribo_cds_counts[["D_S1"]]),
    reads_rep1 = as.numeric(ribo_cds_counts[["D_S1"]]),
    reads_rep2 = as.numeric(ribo_cds_counts[["G_S1"]]),
    mean_read_length = rowMeans(cbind(
      rep(ribo_read_length_means[["D_S1"]], length(ribo_cds_counts[["D_S1"]])),
      rep(ribo_read_length_means[["G_S1"]], length(ribo_cds_counts[["G_S1"]]))
    )),
    condition = "Control"
  )

  # Cap7: F_S3 vs I_S3
  cap7_comparison <- data.frame(
    gene = names(ribo_cds_counts[["F_S3"]]),
    reads_rep1 = as.numeric(ribo_cds_counts[["F_S3"]]),
    reads_rep2 = as.numeric(ribo_cds_counts[["I_S3"]]),
    mean_read_length = rowMeans(cbind(
      rep(ribo_read_length_means[["F_S3"]], length(ribo_cds_counts[["F_S3"]])),
      rep(ribo_read_length_means[["I_S3"]], length(ribo_cds_counts[["I_S3"]]))
    )),
    condition = "Cap7"
  )

  # Combine
  ribo_cds_comparison <- rbind(control_comparison, cap7_comparison)
  ribo_cds_comparison$condition <- factor(ribo_cds_comparison$condition, levels = c("Control", "Cap7"))

  # Calculate correlations
  control_corr <- cor(control_comparison$reads_rep1, control_comparison$reads_rep2, method = "pearson")
  cap7_corr <- cor(cap7_comparison$reads_rep1, cap7_comparison$reads_rep2, method = "pearson")

  cat("[04_qc_correlations] Ribo-seq CDS read correlations:\n")
  cat(sprintf("  Control (D_S1 vs G_S1): r = %.4f\n", control_corr))
  cat(sprintf("  Cap7 (F_S3 vs I_S3):    r = %.4f\n", cap7_corr))

  if (!require("ggplot2", quietly = TRUE)) {
    cat("Warning: ggplot2 not found, skipping Ribo-seq CDS reproducibility plot\n")
  } else {
    p_ribo_cds_repro <- ggplot(ribo_cds_comparison, aes(x = reads_rep1, y = reads_rep2,
                                                         color = mean_read_length, size = reads_rep1 + reads_rep2)) +
      geom_point(alpha = 0.6) +
      facet_wrap(~condition, nrow = 1, scales = "free") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
      scale_color_gradient(low = "#0D0887FF", high = "#E16462FF", name = "Mean read\nlength (nt)") +
      scale_size_continuous(name = "Total CDS\nreads", range = c(1, 4), guide = "none") +
      scale_x_log10() +
      scale_y_log10() +
      theme_classic(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", size = 20),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 14, face = "bold")
      ) +
      labs(
        title = "Ribo-seq CDS reads: reproducibility despite read length heterogeneity",
        subtitle = "High replicate correlation (r > 0.99) validates gene-level analysis despite partial digestion",
        x = "CDS reads - Replicate 1 (log10)",
        y = "CDS reads - Replicate 2 (log10)"
      )

    # Add correlation labels
    p_ribo_cds_repro <- p_ribo_cds_repro +
      annotate("text", x = Inf, y = -Inf,
               label = paste0("Control: r = ", round(control_corr, 4)),
               hjust = 1.1, vjust = -0.5, size = 5, fontface = "bold",
               color = "#0D0887FF") +
      annotate("text", x = Inf, y = -Inf,
               label = paste0("Cap7: r = ", round(cap7_corr, 4)),
               hjust = 1.1, vjust = -1.8, size = 5, fontface = "bold",
               color = "#E16462FF")

    # Export as PNG
    ggsave(
      file.path(dir_results, "ribo_cds_reads_reproducibility.png"),
      plot = p_ribo_cds_repro,
      width = 13, height = 6, dpi = 300
    )
    cat("[04_qc_correlations] Exported PNG: ribo_cds_reads_reproducibility.png\n")

    # Export as SVG
    ggsave(
      file.path(dir_results, "ribo_cds_reads_reproducibility.svg"),
      plot = p_ribo_cds_repro,
      width = 13, height = 6
    )
    cat("[04_qc_correlations] Exported SVG: ribo_cds_reads_reproducibility.svg\n")
  }
}

# ============================================================
# CDS Reads Reproducibility, single panel colored by condition
#
# Same replicate 1 vs replicate 2 comparison as above, but both
# conditions are overlaid in one (unfaceted) panel and
# distinguished only by point color, one plot per assay.
# ============================================================

cat("[04_qc_correlations] Building single-panel CDS reproducibility plots (colored by condition)...\n")

# Display label for the treated condition in these plots (the underlying
# factor level used everywhere else in the pipeline stays "Cap7"; this is
# cosmetic relabeling for these figures only).
cap_label <- "Capreomycin"
condition_colors <- setNames(c("#0D0887FF", "#E16462FF"), c("Control", cap_label))

# Helper: build a rep1-vs-rep2 dataframe from a named-vector counts list
build_replicate_df <- function(counts_list, control_reps, cap_reps, assay_name = NULL) {
  control_df <- data.frame(
    gene = names(counts_list[[control_reps[1]]]),
    reads_rep1 = as.numeric(counts_list[[control_reps[1]]]),
    reads_rep2 = as.numeric(counts_list[[control_reps[2]]]),
    condition = "Control"
  )
  cap_df <- data.frame(
    gene = names(counts_list[[cap_reps[1]]]),
    reads_rep1 = as.numeric(counts_list[[cap_reps[1]]]),
    reads_rep2 = as.numeric(counts_list[[cap_reps[2]]]),
    condition = cap_label
  )
  out <- rbind(control_df, cap_df)
  out$condition <- factor(out$condition, levels = c("Control", cap_label))
  if (!is.null(assay_name)) out$assay <- assay_name
  out
}

plot_replicate_scatter <- function(df, title) {
  p <- ggplot(df, aes(x = reads_rep1, y = reads_rep2, color = condition)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey50", linewidth = 0.8) +
    scale_color_manual(values = condition_colors) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    ) +
    labs(
      title = title,
      x = "CDS reads - Replicate 1 (log10)",
      y = "CDS reads - Replicate 2 (log10)",
      color = "Condition"
    )

  p
}

if (!require("ggplot2", quietly = TRUE)) {
  cat("Warning: ggplot2 not found, skipping single-panel CDS reproducibility plots\n")
} else {
  # RNA-seq
  rna_replicate_df <- build_replicate_df(rna_counts, rna_control_reps, rna_cap_reps, "RNA-seq")
  p_rna_repro <- plot_replicate_scatter(rna_replicate_df, "RNA-seq CDS reads reproducibility")

  ggsave(
    file.path(dir_results, "rna_cds_reads_reproducibility_by_condition.png"),
    plot = p_rna_repro, width = 8, height = 7, dpi = 300
  )
  ggsave(
    file.path(dir_results, "rna_cds_reads_reproducibility_by_condition.svg"),
    plot = p_rna_repro, width = 8, height = 7
  )
  cat("[04_qc_correlations] Exported RNA-seq single-panel reproducibility plot (PNG/SVG)\n")

  # Ribo-seq
  ribo_replicate_df <- build_replicate_df(ribo_counts, ribo_control_reps, ribo_cap_reps, "Ribo-seq")
  p_ribo_repro <- plot_replicate_scatter(ribo_replicate_df, "Ribo-seq CDS reads reproducibility")

  ggsave(
    file.path(dir_results, "ribo_cds_reads_reproducibility_by_condition.png"),
    plot = p_ribo_repro, width = 8, height = 7, dpi = 300
  )
  ggsave(
    file.path(dir_results, "ribo_cds_reads_reproducibility_by_condition.svg"),
    plot = p_ribo_repro, width = 8, height = 7
  )
  cat("[04_qc_correlations] Exported Ribo-seq single-panel reproducibility plot (PNG/SVG)\n")

  # ============================================================
  # Combined figure: RNA-seq and Ribo-seq side by side
  # ============================================================

  combined_replicate_df <- rbind(rna_replicate_df, ribo_replicate_df)
  combined_replicate_df$assay <- factor(combined_replicate_df$assay,
                                         levels = c("RNA-seq", "Ribo-seq"))

  p_combined_repro <- ggplot(combined_replicate_df,
                              aes(x = reads_rep1, y = reads_rep2, color = condition)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey50", linewidth = 0.8) +
    facet_wrap(~assay, nrow = 1, scales = "free") +
    scale_color_manual(values = condition_colors) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 16, face = "bold")
    ) +
    labs(
      title = "CDS reads reproducibility",
      x = "CDS reads - Replicate 1 (log10)",
      y = "CDS reads - Replicate 2 (log10)",
      color = "Condition"
    )

  ggsave(
    file.path(dir_results, "cds_reads_reproducibility_combined.png"),
    plot = p_combined_repro, width = 14, height = 7, dpi = 300
  )
  ggsave(
    file.path(dir_results, "cds_reads_reproducibility_combined.svg"),
    plot = p_combined_repro, width = 14, height = 7
  )
  cat("[04_qc_correlations] Exported combined RNA-seq + Ribo-seq reproducibility plot (PNG/SVG)\n")
}

cat("[04_qc_correlations] Single-panel CDS reproducibility plots complete.\n")
