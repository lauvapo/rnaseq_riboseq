# Conda Environment Setup

Use Conda to manage all dependencies in a single isolated environment. This eliminates version conflicts and ensures reproducibility.

## What is Conda?

Conda is a package and environment manager that:
- ✅ Manages R, Python, and bioinformatics tools in one place
- ✅ Creates isolated environments (no conflicts with system packages)
- ✅ Ensures reproducibility across different machines
- ✅ Simplifies installation of complex bioinformatics tools

## Install Conda

### Option 1: Miniconda (Recommended - minimal)
Lightweight version with just conda, pip, and Python.

```bash
# Download (Linux/macOS)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install
bash Miniconda3-latest-Linux-x86_64.sh

# Restart terminal or run:
source ~/.bashrc
```

### Option 2: Anaconda
Full-featured version (larger download).

Download from: https://www.anaconda.com/download

### Verify Installation

```bash
conda --version
conda list
```

## Create the Capreomycin Environment

Navigate to the project directory and create the environment:

```bash
cd /path/to/paper_capreomycin_remapping

# Create environment from file
conda env create -f environment.yml

# This takes 10-15 minutes (downloads ~3 GB of packages)
```

## Activate the Environment

Before running any scripts, activate the environment:

```bash
conda activate capreomycin-pipeline
```

You should see `(capreomycin-pipeline)` in your terminal prompt.

## Run Pipeline with Conda

```bash
# 1. Activate environment
conda activate capreomycin-pipeline

# 2. Run pipeline
cd /path/to/paper_capreomycin_remapping
Rscript pipeline/run_pipeline.R

# 3. Deactivate when done (optional)
conda deactivate
```

## Included in This Environment

### R (v4.3) + Essential Packages
- **Data manipulation**: dplyr, tidyr, readr, data.table
- **Visualization**: ggplot2, ggrepel, ComplexHeatmap, ggridges, viridis
- **Statistical analysis**: DESeq2 (bioinformatics)
- **Enrichment analysis**: clusterProfiler

### Bioinformatics Tools (for raw data processing)
- **bowtie2** — Read alignment
- **samtools** — BAM/SAM file handling
- **cutadapt** — Adapter trimming
- **subread** — Feature counting (featureCounts)
- **fastqc** — Quality control
- **multiqc** — QC aggregation

### Python (v3.11) + Utilities
- **pandas, numpy, scipy** — Data analysis
- **matplotlib, seaborn** — Plotting
- **biopython** — Sequence utilities
- **pysam** — SAM/BAM file handling

## Common Commands

### Check what's installed
```bash
conda activate capreomycin-pipeline
conda list
```

### Update packages
```bash
conda activate capreomycin-pipeline
conda update --all
```

### Install additional package
```bash
conda activate capreomycin-pipeline
conda install <package_name>
```

### Remove environment
```bash
conda env remove --name capreomycin-pipeline
```

### List all environments
```bash
conda env list
```

## Troubleshooting

### "conda: command not found"
- Miniconda/Anaconda not installed or not in PATH
- Run: `source ~/.bashrc` or restart terminal

### "conda activate doesn't work"
- You need to initialize conda first:
  ```bash
  conda init bash
  # Restart terminal
  ```

### Package conflicts during creation
- Clear conda cache: `conda clean --all`
- Try again: `conda env create -f environment.yml`

### Need specific package version?
Edit `environment.yml` to specify version:
```yaml
- bowtie2=2.5    # Specific version
- samtools>=1.15 # Minimum version
- python=3.11    # Exact version
```

Then update environment:
```bash
conda env update -f environment.yml
```

## Performance Tips

### Speed up package solving
```bash
conda config --set channel_priority strict
```

### Use mamba (faster alternative to conda)
```bash
conda install -c conda-forge mamba
mamba env create -f environment.yml
```

### Avoid disk bloat
```bash
conda clean --all  # Remove unused packages
```

## Share Environment with Others

### Export environment
```bash
conda env export > environment-export.yml
```

### Others can recreate it
```bash
conda env create -f environment-export.yml
```

## Integration with Preprocessing Scripts

The preprocessing scripts (`process_fastq_rnaseq.sh`, `process_fastq_riboseq.sh`) require:
- bowtie2 ✅ (in environment)
- samtools ✅ (in environment)
- cutadapt ✅ (in environment)
- featureCounts (subread) ✅ (in environment)

To use preprocessing scripts with Conda:

```bash
# Activate environment
conda activate capreomycin-pipeline

# Run preprocessing
bash process_fastq_rnaseq.sh
bash process_fastq_riboseq.sh

# This generates counts that feed into the main pipeline
```

## Why Use Conda for This Project?

✅ **Reproducibility**: Same environment on any machine  
✅ **Simplicity**: Install everything with one command  
✅ **Isolation**: No conflicts with system packages  
✅ **Easy sharing**: Send `environment.yml` to collaborators  
✅ **Version control**: Track exact versions used  
✅ **Bioinformatics focus**: Bioconda has thousands of tools  

---

**Next**: Create environment and run pipeline:
```bash
conda env create -f environment.yml
conda activate capreomycin-pipeline
Rscript pipeline/run_pipeline.R
```

**Questions?** Check [Conda documentation](https://docs.conda.io/)
