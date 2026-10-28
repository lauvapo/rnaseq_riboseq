# Capreomycin Effect on M. tuberculosis: RNA-seq and Ribo-seq Analysis Pipeline

A reproducible bioinformatics pipeline for analyzing differential expression and translational response to capreomycin treatment in *Mycobacterium tuberculosis* H37Rv.

## Overview

This pipeline performs comprehensive transcriptomic and translatome analysis:

- **Quality Control**: Replicate correlations, read length distributions
- **Differential Expression**: DESeq2-based analysis for RNA-seq and Ribo-seq
- **Enrichment Analysis**: 
  - Gene Ontology (GO) terms
  - KEGG Pathways
  - Curated regulons
  - iModulons (M. tuberculosis systems biology modules)
- **Visualization**: Volcano plots, heatmaps, correlation matrices, enrichment dotplots

## Project Structure

```
paper_capreomycin_remapping/
├── pipeline/               # Main analysis pipeline
│   ├── 00_config.R        # Configuration & shared parameters
│   ├── run_pipeline.R     # Master driver script
│   ├── 02_quality_control_correlations.R
│   ├── 03_differential_expression.R
│   ├── 04_upset_plot.R
│   ├── 05b_volcano_plots.R
│   ├── 05c_scatter_transcriptome_translatome.R
│   ├── 05d_gene_category_heatmaps.R
│   ├── 05_enrichment_go_pathway.R
│   ├── 05_enrichment_regulon.R
│   ├── 05_enrichment_imodulons.R
│   ├── 05_transcription_network_analysis.R
│   ├── 06_export_spreadsheets.R
│   ├── 07_summary_figures.R
│   ├── 08_enrichment_heatmaps.R
│   ├── 09_enrichment_logfc_heatmaps.R
│   ├── 10_gene_logfc_heatmaps.R
│   ├── 11_curated_gene_heatmap.R
│   └── 12_pathway_ridgeplots.R
├── pathways/               # Annotation files for enrichment analysis
│   ├── gene_map.tsv       # Gene symbol ↔ locus tag mapping
│   ├── go_terms_genes.txt # GO term annotations
│   ├── pathway_genes.txt  # KEGG pathway annotations
│   ├── selected_regulon.gmt # Curated regulon definitions
│   ├── transcrnetwork_modulome.csv # M. tuberculosis modulome network
│   └── imodulons.tsv      # iModulon gene mappings
├── data/                   # Input data (not in GitHub, keep locally)
│   ├── RNA_seq/counts/CDS/ # RNA-seq read counts per gene
│   └── Ribo_seq/counts/CDS/# Ribo-seq read counts per gene
├── results/                # Output files (not in GitHub)
│   ├── DE_and_enrichment_results.xlsx # Excel workbook with all results
│   ├── significative_genes/  # Per-gene DE results
│   ├── enrichment_results/   # ORA tables and visualizations
│   └── *.png, *.svg       # Publication-ready figures
├── scripts/                # Utility scripts
├── .gitignore             # Git ignore rules
├── README.md              # This file
└── INSTALLATION.md        # Setup and installation instructions
```

## Quick Start

### Choose Your Setup Path:

#### 🖥️ **Local Machine with Internet** → Use Conda
```bash
git clone https://github.com/yourusername/capreomycin-pipeline.git
cd capreomycin-pipeline

conda env create -f environment.yml
conda activate capreomycin-pipeline

# Preprocess (optional)
bash pipeline/preprocessing/process_fastq_rnaseq.sh
bash pipeline/preprocessing/process_fastq_riboseq.sh

# Analyze
Rscript pipeline/run_pipeline.R
```
See: **CONDA_SETUP.md** for detailed Conda instructions.

#### 🔧 **Research Server (with/without internet)** → See SERVER_SETUP.md
```bash
git clone https://github.com/yourusername/capreomycin-pipeline.git
cd capreomycin-pipeline

# Automatic installation of R packages
Rscript install_dependencies.R

# Run analysis
Rscript pipeline/run_pipeline.R
```
See: **SERVER_SETUP.md** for server-specific instructions (modules, SLURM, offline setup, etc.)

#### ⚡ **Quick Install** (if you already have count tables)
```bash
git clone https://github.com/yourusername/capreomycin-pipeline.git
cd capreomycin-pipeline

Rscript install_dependencies.R  # Install R packages
Rscript pipeline/run_pipeline.R # Run analysis
```

### Prerequisites by Setup Type

**Local machine with Conda:**
- Conda (Miniconda or Anaconda)
- `environment.yml` handles everything

**Research server:**
- R (≥ 4.0)
- Bioinformatics tools (bowtie2, samtools, cutadapt, featureCounts) - see SERVER_SETUP.md
- `install_dependencies.R` installs R packages
- Raw FASTQ or count tables

### Output

This generates:
- `results/DE_and_enrichment_results.xlsx` — comprehensive results workbook
- PNG and SVG figures in `results/`
- TSV files with detailed statistics

## Input Data Format

### Count Files
Place featureCounts output in:
```
data/RNA_seq/counts/CDS/J_S4_CDS.dedup.csv
data/RNA_seq/counts/CDS/M_S4_CDS.dedup.csv
data/RNA_seq/counts/CDS/L_S6_CDS.dedup.csv
data/RNA_seq/counts/CDS/O_S6_CDS.dedup.csv
data/Ribo_seq/counts/CDS/D_S1_CDS.dedup.csv
data/Ribo_seq/counts/CDS/G_S1_CDS.dedup.csv
data/Ribo_seq/counts/CDS/F_S3_CDS.dedup.csv
data/Ribo_seq/counts/CDS/I_S3_CDS.dedup.csv
```

Each file: gene ID (first column) + counts (subsequent column named after BAM file)

### Annotation Files
All in `pathways/`:
- **gene_map.tsv**: 3-column TSV (name, locus_tag, rv_id)
- **go_terms_genes.txt**: GO annotations (Mycobrowser format)
- **pathway_genes.txt**: KEGG pathway annotations
- **.gmt files**: Gene set definitions for custom enrichment
- **.csv files**: Network and iModulon definitions

## Pipeline Steps

1. **Quality Control** → replicate correlations, read length distributions
2. **Differential Expression** → DESeq2 analysis (RNA-seq + Ribo-seq)
3. **Data Visualization** → volcano plots, scatter plots, heatmaps
4. **Transcription Network Analysis** → modulome topology statistics
5. **Enrichment Analysis** → ORA for GO terms, pathways, regulons, iModulons
6. **Export** → consolidated Excel workbook
7. **Summary Figures** → publication-ready visualizations

## Configuration

Edit `pipeline/00_config.R` to customize:
- Sample names and conditions
- Statistical thresholds (FDR cutoff, logFC cutoff)
- Output directories
- Annotation file paths

## Output Files

### Main Results
- **DE_and_enrichment_results.xlsx** — all results in one workbook:
  - All_DE: full differential expression table
  - DE_summary: count statistics
  - GO sheets: enrichment by BP/CC/MF/Pathway
  - Regulon sheets: regulon-based enrichment
  - iModulon sheets: iModulon enrichment

### Visualizations
- **volcano_combined.png/svg** — RNA-seq vs Ribo-seq volcano plots
- **heatmap_*.png/svg** — Gene category heatmaps (lipid, ribosome, translation, transport, etc.)
- ***_enrichment_heatmap.png/svg** — Top enriched terms per gene set
- **ridgeplot_*.png/svg** — Pathway distribution ridge plots
- **summary_*.png** — Publication-ready summary figures

### Detailed Results
- `significative_genes/all_info.tsv` — all DE genes with statistics
- `enrichment_results/*_ora_results.tsv` — ORA tables for each enrichment type
- `enrichment_results/network/` — network topology statistics

## Methods

### Differential Expression
- **Tool**: DESeq2
- **Design**: Two-factor design (Condition: Control vs. Capreomycin)
- **Significance**: p.adjust < 0.1 AND |logFC| > 0.5
- **Output**: normalized counts, logFC, p-values, adjusted p-values

### Over-Representation Analysis (ORA)
- **Tool**: clusterProfiler::enricher()
- **Test**: Hypergeometric test
- **Multiple testing correction**: Benjamini-Hochberg (FDR)
- **Significance**: p.adjust < 0.1, gene set size ≥ 3

### Gene Sets
1. **GO terms**: 3 aspects (Biological Process, Cellular Component, Molecular Function)
2. **Pathways**: KEGG pathway definitions
3. **Regulons**: Curated transcription factor regulons (M. tuberculosis)
4. **iModulons**: Systems biology modules from M. tuberculosis modulome

## Citation

If you use this pipeline, please cite:

```bibtex
[Citation information to be updated]
```

## Author & Contact

**Developer**: [Your name]  
**Institution**: [Your institution]  
**Email**: [Your email]

## License

[Specify license: MIT, GPL-3.0, CC-BY-4.0, etc.]

## Troubleshooting

### Missing input files
- Check `data/RNA_seq/counts/CDS/` and `data/Ribo_seq/counts/CDS/` exist
- Verify sample names match those in `pipeline/00_config.R`
- Ensure annotation files are in `pathways/`

### R package errors
- Run `install.packages()` for missing packages (see INSTALLATION.md)
- Check R version ≥ 4.0

### Low or no enrichment results
- Verify FDR and logFC cutoffs in `pipeline/00_config.R`
- Check gene set sizes (minimum 3 genes required)
- Confirm annotation file format matches expectations

### Out of memory
- Reduce figure resolution (dpi parameter in ggsave calls)
- Process in smaller batches if available

## Related Files

- **INSTALLATION.md** — Detailed setup instructions
- **pipeline/README.md** — Technical documentation for each pipeline step
- **pathways/README.md** — Annotation file documentation

---

**Last updated**: 2026-10-28  
**Pipeline version**: 1.0
