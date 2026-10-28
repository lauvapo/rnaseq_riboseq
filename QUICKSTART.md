# Quick Start Guide

Get the pipeline running in 5 minutes.

## Prerequisites

- R 4.0+
- Raw count files in `data/RNA_seq/counts/CDS/` and `data/Ribo_seq/counts/CDS/`
- R packages installed (see INSTALLATION.md if needed)

## Run in 3 Steps

### 1. Prepare Your Data

Place your featureCounts output files in the correct directories:

```
data/RNA_seq/counts/CDS/
├── J_S4_CDS.dedup.csv
├── M_S4_CDS.dedup.csv
├── L_S6_CDS.dedup.csv
└── O_S6_CDS.dedup.csv

data/Ribo_seq/counts/CDS/
├── D_S1_CDS.dedup.csv
├── G_S1_CDS.dedup.csv
├── F_S3_CDS.dedup.csv
└── I_S3_CDS.dedup.csv
```

### 2. Check Configuration (Optional)

Edit `pipeline/00_config.R` if needed to change:
- Sample names
- FDR or logFC thresholds
- Output directory

### 3. Run the Pipeline

```bash
Rscript pipeline/run_pipeline.R
```

**Done!** Results appear in `results/` after 30-60 minutes.

## Check Results

- **Comprehensive results**: `results/DE_and_enrichment_results.xlsx`
- **Figures**: `results/*.png` and `results/*.svg`
- **Detailed tables**: `results/enrichment_results/`

## Common Next Steps

**Adjust thresholds?**
```r
# In pipeline/00_config.R:
fdr_cutoff <- 0.05    # Stricter
logfc_cutoff <- 1.0
```

**Run only certain steps?**
Edit `pipeline/run_pipeline.R` and comment out unwanted steps.

**Change sample names?**
Edit the sample vectors in `pipeline/00_config.R`:
```r
rna_samples <- c("your_sample1", "your_sample2", ...)
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "File not found" | Check `data/` directory paths in 00_config.R |
| "Package X not found" | Install missing package: `install.packages("X")` |
| "No enrichment results" | Reduce fdr_cutoff or logfc_cutoff thresholds |
| "Out of memory" | Use a machine with more RAM or close other apps |

## Full Documentation

- **README.md** — Complete project documentation
- **INSTALLATION.md** — Detailed setup instructions
- **pipeline/README.md** — Technical details for each step

---

Ready to run? Just execute:
```bash
Rscript pipeline/run_pipeline.R
```
