# ============================================================
# Step 5a: Transcription Network Analysis
#
# Analyzes the transcription network from transcrnetwork_modulome.csv
# to identify key regulators, target distributions, and network
# topology. Independent of enrichment analyses.
#
# Outputs:
#   enrichment_results/network/network_regulators.tsv
#   enrichment_results/network/network_targets.tsv
#   enrichment_results/network/network_topology.tsv
#   enrichment_results/network/ (PNG visualizations)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(igraph)
})

# Load transcription network (Modulome format)
transcrnetwork <- read.csv(transcrnetwork_file, stringsAsFactors = FALSE)
cat("[05a_network] Loaded transcription network with", nrow(transcrnetwork),
    "regulator-target edges\n")

# Create output directory
out_dir <- file.path(dir_enrichment, "network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. REGULATOR ANALYSIS
# ============================================================
# Calculate statistics for each regulator
regulator_stats <- transcrnetwork %>%
  group_by(regulator_id) %>%
  summarize(
    num_targets = n_distinct(target_id),
    source_count = n_distinct(source),
    .groups = "drop"
  ) %>%
  arrange(desc(num_targets))

write.table(regulator_stats, file.path(out_dir, "network_regulators.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[05a_network] Regulator statistics:\n")
cat("[05a_network]   Total regulators:", nrow(regulator_stats), "\n")
cat("[05a_network]   Mean targets per regulator:",
    round(mean(regulator_stats$num_targets), 2), "\n")
cat("[05a_network]   Median targets per regulator:",
    median(regulator_stats$num_targets), "\n")
cat("[05a_network]   Top 5 regulators by target count:\n")
top_regs <- regulator_stats %>% head(5)
for (i in 1:nrow(top_regs)) {
  cat("[05a_network]    ", top_regs$regulator_id[i], "-",
      top_regs$num_targets[i], "targets\n")
}

# ============================================================
# 2. TARGET ANALYSIS
# ============================================================
# Calculate statistics for each target gene
target_stats <- transcrnetwork %>%
  group_by(target_id) %>%
  summarize(
    num_regulators = n_distinct(regulator_id),
    .groups = "drop"
  ) %>%
  arrange(desc(num_regulators))

write.table(target_stats, file.path(out_dir, "network_targets.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[05a_network] Target gene statistics:\n")
cat("[05a_network]   Total target genes:", nrow(target_stats), "\n")
cat("[05a_network]   Mean regulators per gene:",
    round(mean(target_stats$num_regulators), 2), "\n")
cat("[05a_network]   Median regulators per gene:",
    median(target_stats$num_regulators), "\n")
cat("[05a_network]   Top 5 genes by regulator count:\n")
top_targets <- target_stats %>% head(5)
for (i in 1:nrow(top_targets)) {
  cat("[05a_network]    ", top_targets$target_id[i], "-",
      top_targets$num_regulators[i], "regulators\n")
}

# ============================================================
# 3. NETWORK TOPOLOGY ANALYSIS
# ============================================================
# Create igraph object for network analysis
edges <- transcrnetwork %>%
  dplyr::select(regulator_id, target_id) %>%
  as.data.frame()

g <- graph_from_data_frame(edges, directed = TRUE)

# Calculate network metrics
topology_metrics <- data.frame(
  metric = c(
    "Total edges",
    "Total regulators",
    "Total targets",
    "Total unique genes",
    "Network density",
    "Average degree",
    "Number of connected components"
  ),
  value = c(
    ecount(g),
    length(unique(transcrnetwork$regulator_id)),
    length(unique(transcrnetwork$target_id)),
    length(V(g)),
    edge_density(g),
    mean(degree(g)),
    components(g)$no
  )
)

write.table(topology_metrics, file.path(out_dir, "network_topology.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("[05a_network] Network topology:\n")
cat("[05a_network]   Total edges:", ecount(g), "\n")
cat("[05a_network]   Network density:", round(edge_density(g), 4), "\n")
cat("[05a_network]   Average degree:", round(mean(degree(g)), 2), "\n")
cat("[05a_network]   Connected components:", components(g)$no, "\n")

# ============================================================
# 4. SOURCE ANALYSIS
# ============================================================
source_stats <- transcrnetwork %>%
  group_by(source) %>%
  summarize(
    count = n(),
    num_unique_regulators = n_distinct(regulator_id),
    num_unique_targets = n_distinct(target_id),
    .groups = "drop"
  )

write.table(source_stats, file.path(out_dir, "network_source_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[05a_network] Edge source distribution:\n")
for (i in 1:nrow(source_stats)) {
  cat("[05a_network]   ", source_stats$source[i], ":",
      source_stats$count[i], "edges\n")
}

cat("[05a_network] Done. Network analysis results written to", out_dir, "\n")
