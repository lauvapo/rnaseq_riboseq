#!/usr/bin/env Rscript
# ============================================================
# Automated dependency installer for servers without internet
# Usage: Rscript install_dependencies.R
# ============================================================

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║  Installing Capreomycin Pipeline R Dependencies          ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# CRAN packages (data manipulation, visualization, utilities)
cran_packages <- c(
  # Visualization
  "ggplot2",
  "ggrepel",
  "gridExtra",
  "ggridges",
  "viridis",
  "scales",
  "cowplot",

  # Data manipulation
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "openxlsx",
  "stringr",
  "forcats",
  "lubridate",
  "glue",
  "data.table",

  # Development & utilities
  "devtools",
  "rmarkdown",
  "tidyverse"
)

# Bioconductor packages
bioc_packages <- c(
  "DESeq2",
  "clusterProfiler",
  "ComplexHeatmap"
)

# Install CRAN packages
cat("📦 Installing CRAN packages...\n")
cat("   Repository: http://cran.r-project.org\n\n")

for (pkg in cran_packages) {
  cat(sprintf("   Installing %-30s", pkg))

  tryCatch({
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg,
                      repos = "http://cran.r-project.org",
                      dependencies = TRUE,
                      verbose = FALSE,
                      quiet = TRUE)
      cat(" ✓\n")
    } else {
      cat(" (already installed) ✓\n")
    }
  }, error = function(e) {
    cat(sprintf(" ✗ FAILED\n   Error: %s\n", e$message))
  })
}

# Install Bioconductor packages
cat("\n📦 Installing Bioconductor packages...\n")
cat("   Repository: Bioconductor\n\n")

tryCatch({
  if (!require("BiocManager", quietly = TRUE)) {
    cat("   Installing BiocManager...\n")
    install.packages("BiocManager",
                    repos = "http://cran.r-project.org",
                    verbose = FALSE,
                    quiet = TRUE)
  }
}, error = function(e) {
  cat(sprintf("   ✗ Failed to install BiocManager: %s\n", e$message))
})

for (pkg in bioc_packages) {
  cat(sprintf("   Installing %-30s", pkg))

  tryCatch({
    if (!require(pkg, character.only = TRUE)) {
      suppressMessages({
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      })
      cat(" ✓\n")
    } else {
      cat(" (already installed) ✓\n")
    }
  }, error = function(e) {
    cat(sprintf(" ✗ FAILED\n   Error: %s\n", e$message))
  })
}

# Verify installations
cat("\n✅ Verifying installations...\n\n")

all_packages <- c(cran_packages, bioc_packages)
installed <- sapply(all_packages, function(pkg) {
  require(pkg, character.only = TRUE, quietly = TRUE)
})

if (all(installed)) {
  cat("✓ All packages installed successfully!\n")
  cat("\n╔════════════════════════════════════════════════════════════╗\n")
  cat("║  ✓ Ready to run pipeline!                                 ║\n")
  cat("║                                                            ║\n")
  cat("║  Next step: Rscript pipeline/run_pipeline.R              ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n")
  quit(status = 0)
} else {
  failed <- names(installed[!installed])
  cat(sprintf("\n✗ Failed to install %d packages:\n", length(failed)))
  for (pkg in failed) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nPlease install manually or contact your administrator.\n")
  quit(status = 1)
}
