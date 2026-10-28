# Server Setup Guide

Instructions for setting up the pipeline on a research server **with or without internet access**.

## Quick Start (If Server Has Internet)

```bash
git clone https://github.com/yourusername/capreomycin-pipeline.git
cd capreomycin-pipeline

conda env create -f environment.yml
conda activate capreomycin-pipeline

Rscript pipeline/run_pipeline.R
```

---

## Setup on Server WITHOUT Internet Access

### Approach 1: Pre-downloaded Environment (Recommended)

**Step 1: Download environment on a connected machine**

```bash
# On any machine with internet and conda installed
git clone https://github.com/yourusername/capreomycin-pipeline.git
cd capreomycin-pipeline

# Download and pack environment
conda env create --download-only -f environment.yml
conda pack -n capreomycin-pipeline -o capreomycin-pipeline.tar.gz

# Compress further for transfer
gzip capreomycin-pipeline.tar.gz
```

**Step 2: Transfer to server**

```bash
# On your local machine
scp capreomycin-pipeline.tar.gz user@server:/path/to/destination/

# Or use rsync (faster for large files)
rsync -avz capreomycin-pipeline.tar.gz user@server:/path/to/destination/
```

**Step 3: Unpack on server**

```bash
# On server
cd /path/to/destination/

# Create conda envs directory if needed
mkdir -p ~/miniconda3/envs

# Extract environment
tar -xzf capreomycin-pipeline.tar.gz -C ~/miniconda3/envs/capreomycin-pipeline

# Activate
conda activate capreomycin-pipeline

# Verify
which Rscript
which bowtie2
```

---

### Approach 2: Manual R Installation (If Conda Unavailable)

If your server doesn't have conda or you can't transfer the environment:

**Step 1: Clone repository**

```bash
git clone https://github.com/yourusername/capreomycin-pipeline.git
cd capreomycin-pipeline
```

**Step 2: Install R packages manually**

Create a file `install_r_packages.R`:

```r
# Install CRAN packages
packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "readxl", "openxlsx",
  "ggrepel", "gridExtra", "ggridges", "viridis", "scales", "cowplot",
  "stringr", "forcats", "lubridate", "glue", "data.table", "devtools"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org")
  }
}

# Install Bioconductor packages
if (!require("BiocManager")) {
  install.packages("BiocManager", repos = "http://cran.r-project.org")
}

bioc_packages <- c("DESeq2", "clusterProfiler", "ComplexHeatmap")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

cat("All R packages installed successfully!\n")
```

Run it:

```bash
Rscript install_r_packages.R
```

**Step 3: Check for bioinformatics tools**

```bash
# Check if these tools are already on server
which bowtie2
which samtools
which cutadapt
which featureCounts

# If missing, ask your server admin to install them, or load them:
module load bowtie2
module load samtools
module load cutadapt
module load subread  # for featureCounts
```

**Step 4: Run pipeline**

```bash
Rscript pipeline/run_pipeline.R
```

---

## Checking Server Environment

### Check What's Available

```bash
# Available modules (if using module system)
module avail

# List installed tools
which R
which Rscript
which bowtie2
which samtools
which cutadapt
which featureCounts

# Check conda
conda --version
conda env list
```

### Load Required Tools (Module System)

If your server uses modules:

```bash
# Load R
module load r/4.3

# Load bioinformatics tools
module load bowtie2/2.5
module load samtools/1.17
module load cutadapt/4.4
module load subread/2.0

# Verify
Rscript --version
bowtie2 --version
samtools --version
```

### Create a `.bashrc` alias (Optional)

Add to your `~/.bashrc` to auto-load dependencies:

```bash
# Add this to ~/.bashrc
activate_pipeline() {
    module load r/4.3
    module load bowtie2/2.5
    module load samtools/1.17
    module load cutadapt/4.4
    module load subread/2.0
    conda activate capreomycin-pipeline
    echo "Pipeline environment activated!"
}
```

Then activate with:

```bash
source ~/.bashrc
activate_pipeline
```

---

## Running on Server with SLURM/Job Scheduler

### Create a job script (`run_pipeline.sh`)

```bash
#!/bin/bash
#SBATCH --job-name=capreomycin-pipeline
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=pipeline_%j.log

# Load modules
module load r/4.3
module load bowtie2/2.5
module load samtools/1.17
module load cutadapt/4.4
module load subread/2.0

# Activate conda environment (if using conda)
conda activate capreomycin-pipeline

# Change to project directory
cd /path/to/capreomycin-pipeline

# Run pipeline
Rscript pipeline/run_pipeline.R
```

### Submit job

```bash
sbatch run_pipeline.sh

# Check status
squeue -u $USER
```

---

## Troubleshooting Server Setup

### Problem: "R: command not found"

**Solution**: Load R module or install it

```bash
module load r/4.3
# or
Rscript --version  # Check if in PATH
```

### Problem: "bowtie2: command not found"

**Solution**: Install or load bowtie2

```bash
module load bowtie2/2.5
# or ask admin to install
```

### Problem: Package not found in R

**Solution**: Install from CRAN

```r
install.packages("package_name", repos = "http://cran.r-project.org")
```

### Problem: Conda environment can't find packages

**Solution**: Re-create environment on connected machine and transfer tar.gz

```bash
# On connected machine
conda pack -n capreomycin-pipeline -o capreomycin-pipeline.tar.gz

# Transfer and unpack on server
scp capreomycin-pipeline.tar.gz user@server:~/
tar -xzf ~/capreomycin-pipeline.tar.gz -C ~/miniconda3/envs/capreomycin-pipeline
```

### Problem: Out of memory during pipeline

**Solution**: Reduce threads or increase available memory

Edit `pipeline/00_config.R`:

```r
# Reduce threads
THREADS = 4  # Instead of 8
```

Or submit longer SLURM job with more memory:

```bash
#SBATCH --mem=64G
#SBATCH --time=8:00:00
```

---

## Server-Specific Configuration

### Edit config for server paths

If your server has specific paths for data:

```r
# In pipeline/00_config.R:
project_dir <- "/data_shared/lvaron/rnaseq_riboseq"  # Server path
dir_rna_counts_cds <- "/data_shared/raw_counts/RNA"  # Server data location
```

### Store results on fast storage

```r
# Store large results on fast filesystem
dir_results <- "/scratch/lvaron/results"  # Fast SSD
```

---

## Performance Optimization for Servers

### Parallel Processing

```r
# In pipeline/00_config.R, enable parallel processing
library(parallel)
num_cores <- detectCores() - 1
```

### Memory-Efficient Mode

For servers with limited memory:

```r
# In differential_expression script
# Use lower count thresholds
min_count <- 1  # Instead of 10
```

---

## Deployment Checklist

- [ ] Clone repo on server
- [ ] Check R version: `Rscript --version`
- [ ] Check bioinformatics tools installed
- [ ] Load necessary modules
- [ ] Install missing R packages: `Rscript install_r_packages.R`
- [ ] Update `pipeline/00_config.R` with server paths
- [ ] Test with small subset of data
- [ ] Run full pipeline with SLURM job

---

## Example: Complete Server Setup Session

```bash
# 1. Clone repo
git clone https://github.com/yourusername/capreomycin-pipeline.git
cd capreomycin-pipeline

# 2. Check environment
module avail r/
module load r/4.3
Rscript --version

# 3. Install R packages
Rscript install_r_packages.R

# 4. Load bioinformatics tools
module load bowtie2/2.5
module load samtools/1.17
module load cutadapt/4.4
module load subread/2.0

# 5. Verify all tools
which bowtie2 samtools cutadapt featureCounts Rscript

# 6. Run pipeline
Rscript pipeline/run_pipeline.R

# Or submit as job
sbatch run_pipeline.sh
```

---

## Contact & Support

If encountering server-specific issues:

1. Check **SERVER_SETUP.md** (this file)
2. Check **INSTALLATION.md** for manual R package installation
3. Ask your server admin about available modules
4. Check pipeline logs: `tail -f pipeline_*.log`

---

**Next**: Set up your specific server environment using the approach above, then run the pipeline!
