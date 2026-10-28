# Installation & Setup Guide

Complete step-by-step instructions to set up and run the capreomycin pipeline.

## System Requirements

- **OS**: Linux, macOS, or Windows (with WSL2)
- **R**: Version 4.0 or later
- **Disk Space**: ~5 GB (raw data) + ~2 GB (results)
- **RAM**: 8 GB minimum (16 GB recommended)

## Step 1: Install Conda (Recommended)

Conda manages R, Python, bioinformatics tools, and all dependencies in one isolated environment.

### Install Miniconda (Lightweight - Recommended)

```bash
# Download (Linux)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install
bash Miniconda3-latest-Linux-x86_64.sh

# Restart terminal
source ~/.bashrc
```

### Or Install Anaconda
Download from: https://www.anaconda.com/download

### Verify Installation
```bash
conda --version
conda list
```

## Step 2: Clone or Download the Repository

```bash
# Clone from GitHub (once repository is available)
git clone https://github.com/yourusername/capreomycin-pipeline.git
cd capreomycin-pipeline

# Or download and extract the ZIP file manually
```

## Step 3: Create Conda Environment (Recommended)

This is the easiest way to install all dependencies (R, bioinformatics tools, Python, and all packages):

```bash
# From project directory
conda env create -f environment.yml

# This takes 10-15 minutes (downloads ~3 GB)
```

Every time you work on this project, activate the environment:

```bash
conda activate capreomycin-pipeline

# You should see (capreomycin-pipeline) in your terminal prompt
```

**That's it for setup!** Skip to **Step 7: Run the Pipeline**.

For more details on Conda, see **CONDA_SETUP.md**.

---

## Alternative: Manual Installation (Not Recommended)

Skip this if you used Conda above.

### Install R Manually

#### Linux (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install r-base r-base-dev
```

#### macOS
```bash
brew install r
```

#### Windows
Download from: https://cran.r-project.org/

---

## Step 4: Prepare Input Data

### Directory Structure
Create and populate the `data/` directory:

```
data/
├── RNA_seq/
│   └── counts/
│       └── CDS/
│           ├── J_S4_CDS.dedup.csv    (Control replicate 1)
│           ├── M_S4_CDS.dedup.csv    (Control replicate 2)
│           ├── L_S6_CDS.dedup.csv    (Capreomycin replicate 1)
│           └── O_S6_CDS.dedup.csv    (Capreomycin replicate 2)
└── Ribo_seq/
    └── counts/
        └── CDS/
            ├── D_S1_CDS.dedup.csv    (Control replicate 1)
            ├── G_S1_CDS.dedup.csv    (Control replicate 2)
            ├── F_S3_CDS.dedup.csv    (Capreomycin replicate 1)
            └── I_S3_CDS.dedup.csv    (Capreomycin replicate 2)
```

### Data Format
Each count file should be a featureCounts output with:
- Column 1: Gene ID (e.g., "Rv0001")
- Column 2: Gene name (e.g., "dnaA")
- Column 3+: Count columns (last column named after BAM file, e.g., "./alignments/J_S4_CDS.dedup.bam")

Example:
```
Geneid    GeneName    Chromosome    Start    End    Strand    Length    ./alignments/J_S4_CDS.dedup.bam
Rv0001    dnaA        Chromosome    164      1472   +         1309      145
Rv0002    dnaN        Chromosome    1766     2797   -         1032      89
```

## Step 4: Verify Annotation Files

The `pathways/` directory includes required annotation files:

```
pathways/
├── gene_map.tsv                    (gene symbols ↔ locus tags mapping)
├── go_terms_genes.txt              (Gene Ontology annotations)
├── pathway_genes.txt               (KEGG pathway annotations)
├── selected_regulon.gmt            (Regulon definitions - GMT format)
├── transcrnetwork_modulome.csv     (M. tuberculosis regulome network)
└── imodulons.tsv                   (iModulon gene mappings)
```

**All annotation files are included in the repository.**

If you need to update them:
- **GO and pathway files**: Download from [Mycobrowser](https://mycobrowser.epfl.ch/)
- **Regulons**: Download from [RegulonDB](http://regulondb.ccg.unam.mx/) or [DOOR](http://csbl.bmb.uga.edu/DOOR/)
- **iModulons**: Available from supplementary data of Arrieta-Ortiz et al. (2020)

## Step 5: Install R Dependencies

### Quick Installation (Recommended)

Run this in R:

```r
# Install required packages
packages <- c(
  "DESeq2",        # Differential expression analysis
  "ggplot2",       # Plotting
  "dplyr",         # Data manipulation
  "tidyr",         # Data tidying
  "readxl",        # Read Excel files
  "openxlsx",      # Write Excel files
  "clusterProfiler", # Enrichment analysis
  "org.Mm.eg.db",  # Gene annotation (M. tuberculosis custom required)
  "ComplexHeatmap", # Advanced heatmaps
  "ggrepel",       # Intelligent text labels
  "gridExtra",     # Grid layouts
  "ggridges",      # Ridge plots
  "viridis"        # Color scales
)

# Install from CRAN
install.packages(packages)

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("ComplexHeatmap")
```

### Or Use the Installation Script

```bash
Rscript scripts/install_dependencies.R
```

### Verify Installation

```r
# In R console, check all packages load without error
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(openxlsx)
# etc.
```

## Step 6: Configure the Pipeline

Edit `pipeline/00_config.R` to customize:

```r
# Project directory (absolute path)
project_dir <- "/home/user/capreomycin-pipeline"

# Sample metadata (if different from defaults)
rna_samples <- c("J_S4", "M_S4", "L_S6", "O_S6")
ribo_samples <- c("D_S1", "G_S1", "F_S3", "I_S3")

# Statistical thresholds
fdr_cutoff <- 0.1      # 10% false discovery rate
logfc_cutoff <- 0.5    # ~1.4x fold-change minimum

# ORA parameters
ora_pvalue_cutoff <- 0.1   # Enrichment significance
ora_min_gs_size <- 3       # Minimum gene set size
```

## Step 7: Run the Pipeline

### Full Pipeline
```bash
cd /path/to/capreomycin-pipeline
Rscript pipeline/run_pipeline.R
```

This runs all 12 analysis steps sequentially (~30-60 minutes depending on system).

### Run Individual Steps (Optional)
```bash
# Source config first
Rscript -e "source('pipeline/00_config.R')"

# Then run individual steps
Rscript -e "source('pipeline/00_config.R'); source('pipeline/02_quality_control_correlations.R')"
Rscript -e "source('pipeline/00_config.R'); source('pipeline/03_differential_expression.R')"
# etc.
```

## Step 8: Check Results

After the pipeline completes:

1. **Excel workbook**: `results/DE_and_enrichment_results.xlsx`
   - Contains all DE results and enrichment tables
   - Multiple sheets: GO, pathways, regulons, iModulons

2. **Figures**: `results/*.png` and `results/*.svg`
   - Volcano plots
   - Heatmaps
   - Enrichment visualizations
   - Summary figures

3. **Detailed results**: `results/enrichment_results/*.tsv`
   - ORA tables with p-values and gene lists

## Troubleshooting

### "Package X not found"
```r
# Install the missing package
install.packages("X")
# Or from Bioconductor
BiocManager::install("X")
```

### "File not found" errors
- Check all paths in `pipeline/00_config.R` are absolute and correct
- Verify input files exist with `ls -la data/RNA_seq/counts/CDS/`
- Check annotation files in `pathways/`

### "No such file" during featureCounts read
- Verify column names in count files match expected format
- Check file encoding (should be UTF-8)
- Ensure no blank rows/columns in count files

### Out of memory during DESeq2
- Close other applications
- Increase system RAM or use a machine with more memory
- Note: Pipeline should handle standard 16 GB RAM

### Enrichment results are empty
- Check FDR/logFC thresholds aren't too stringent
- Verify gene symbols in annotation files match gene_map.tsv
- Confirm GO/pathway annotation files aren't malformed

### Plotting errors
- Ensure `results/` directory is writable
- Check disk space available
- Verify all required fonts are installed

## Advanced Configuration

### Custom Thresholds

Edit `pipeline/00_config.R`:
```r
# Stricter thresholds (fewer significant genes)
fdr_cutoff <- 0.05
logfc_cutoff <- 1.0

# Looser thresholds (more significant genes)
fdr_cutoff <- 0.2
logfc_cutoff <- 0.3
```

### Custom Gene Sets

Add new .gmt files to `pathways/` and reference in pipeline scripts:
```r
custom_geneset <- read.delim("pathways/my_custom_set.gmt", ...)
```

### Output Directory

Change where results are saved:
```r
dir_results <- "/path/to/custom/results"
```

## Performance Optimization

### Speed Up Pipeline
- **Parallelize**: Modify pipeline scripts to use `parallel::mclapply()`
- **Skip steps**: Comment out unnecessary steps in `run_pipeline.R`
- **Reduce resolution**: Lower dpi in ggsave() calls

### Reduce Memory Usage
- Process samples in batches
- Use sparse matrix representations for large datasets
- Export intermediate results to free RAM

## Getting Help

1. **Check logs**: Pipeline prints detailed messages to console
2. **Review comments**: Each pipeline script has inline documentation
3. **Verify data**: Ensure input files match expected format
4. **Check GitHub Issues**: See if problem is already documented
5. **Contact author**: Email for specific issues

## Next Steps

Once installation is complete:
1. Read **README.md** for project overview
2. Review **pipeline/README.md** for technical details
3. Run the full pipeline: `Rscript pipeline/run_pipeline.R`
4. Explore results in `results/DE_and_enrichment_results.xlsx`

---

**Installation version**: 1.0  
**Last updated**: 2026-10-28
