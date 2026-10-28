# Pipeline: RNA-seq & Ribo-seq Differential Expression & Enrichment Analysis

Automated end-to-end analysis pipeline for capreomycin-induced transcriptomic and translatome changes in *M. tuberculosis* H37Rv.

## Overview

This pipeline processes raw featureCounts output through:
1. Quality control (replicate correlations)
2. Differential expression analysis (DESeq2)
3. Data visualization (volcano plots, heatmaps, scatter plots)
4. Transcription network analysis (modulome topology)
5. Enrichment analysis (GO, pathways, regulons, iModulons)
6. Export to Excel workbook and publication-ready figures

## Prerequisites

**Input Data**: Raw count files from featureCounts in:
- `../data/RNA_seq/counts/CDS/*_CDS.dedup.csv` (RNA-seq)
- `../data/Ribo_seq/counts/CDS/*_CDS.dedup.csv` (Ribo-seq)

See **../INSTALLATION.md** for data preparation.

## Run the Pipeline

```bash
cd /path/to/paper_capreomycin_remapping
Rscript pipeline/run_pipeline.R
```

Pipeline takes 30-60 minutes depending on system resources.

## Pipeline Steps

| # | Script | Description | Key Outputs |
|---|--------|-------------|------------|
| 00 | `00_config.R` | Load configuration, paths, sample metadata, annotation files | (sourced by all steps) |
| 02 | `02_quality_control_correlations.R` | Replicate correlations, RNA-seq vs Ribo-seq correlation, read length distributions | `results/replicate_correlation.tsv`, `results/read_length_summary.tsv` |
| 03 | `03_differential_expression.R` | DESeq2 analysis for RNA-seq and Ribo-seq; classify genes by regulation pattern | `significative_genes/all_info.tsv`, `significative_genes/*_sig.tsv` |
| 04 | `04_upset_plot.R` | UpSet plot showing overlap of up/down-regulated genes across assays | `results/upset_rna_ribo.png`, `.svg` |
| 05b | `05b_volcano_plots.R` | Volcano plots (separate + combined) with labeled key genes | `results/volcano_combined.png`, `volcano_rna.png`, `volcano_ribo.png` |
| 05c | `05c_scatter_transcriptome_translatome.R` | Scatter plot: RNA-seq logFC vs Ribo-seq logFC | `results/scatter_translatome.png`, `.svg` |
| 05d | `05d_gene_category_heatmaps.R` | Curated heatmaps by biological theme (translation, lipids, transport, etc.) | `results/curated_theme_*.png`, `.svg` |
| 05 | `05_transcription_network_analysis.R` | M. tuberculosis modulome network topology analysis | `enrichment_results/network/*.tsv` |
| 05 | `05_enrichment_go_pathway.R` | ORA for GO terms and KEGG pathways | `enrichment_results/go_pathway_ora_results.tsv` |
| 05 | `05_enrichment_regulon.R` | ORA for curated regulons | `enrichment_results/regulon_ora_results.tsv` |
| 05 | `05_enrichment_imodulons.R` | ORA for M. tuberculosis iModulons | `enrichment_results/imodulons_ora_results.tsv` |
| 06 | `06_export_spreadsheets.R` | Consolidate all results into single Excel workbook | `results/DE_and_enrichment_results.xlsx` |
| 07 | `07_summary_figures.R` | Publication-ready summary figures and captions | `results/summary_*.png`, `results/figure_captions.txt` |
| 08 | `08_enrichment_heatmaps.R` | Term × gene-set ORA heatmaps | `enrichment_results/*_heatmap.png`, `.svg` |
| 09 | `09_enrichment_logfc_heatmaps.R` | Mean logFC per term for significant ORA hits | `enrichment_results/logfc_heatmap*.png`, `.svg` |
| 10 | `10_gene_logfc_heatmaps.R` | Gene-level logFC heatmap by DE classification | `results/gene_logfc_heatmap.png`, `.svg` |
| 11 | `11_curated_gene_heatmap.R` | Curated gene heatmap grouped by biological theme | `results/curated_gene_heatmap.png`, `.svg` |
| 12 | `12_pathway_ridgeplots.R` | Ridge plots of top pathway categories (separate + combined) | `results/ridgeplot_*.png`, `.svg` |

## Configuration

All parameters are in **`00_config.R`**:

```r
# Sample names and conditions
rna_samples  <- c("J_S4", "M_S4", "L_S6", "O_S6")
ribo_samples <- c("D_S1", "G_S1", "F_S3", "I_S3")

# Statistical thresholds
fdr_cutoff   <- 0.1      # p.adjust < 0.1
logfc_cutoff <- 0.5      # |logFC| > 0.5

# Over-representation analysis parameters
ora_pvalue_cutoff <- 0.1
ora_qvalue_cutoff <- 0.1
ora_min_gs_size   <- 3   # Minimum genes per set
```

To customize:
1. Edit `00_config.R`
2. Re-run `Rscript pipeline/run_pipeline.R`

## Output Files

### Main Results Workbook
- **`results/DE_and_enrichment_results.xlsx`** — All results in one workbook:
  - **All_DE**: Full differential expression table with statistics
  - **DE_summary**: Counts per condition/assay
  - **Replicate_correlation**: Sample correlation matrix
  - **GO_BP_*, GO_CC_*, GO_MF_***: Gene Ontology enrichment
  - **Pathway_***: KEGG pathway enrichment
  - **Regulon_***: Regulon-based enrichment
  - **iModulon_***: iModulon enrichment

### Detailed Results
```
significative_genes/
├── all_info.tsv          # All genes with DE statistics
├── rna_sig.tsv           # Significant RNA-seq genes
├── ribo_sig.tsv          # Significant Ribo-seq genes
└── *_gene_counts.tsv     # Category breakdowns

enrichment_results/
├── go_pathway_ora_results.tsv    # GO/pathway ORA table
├── regulon_ora_results.tsv       # Regulon ORA table
├── imodulons_ora_results.tsv     # iModulon ORA table
├── network/                      # Network topology analysis
│   ├── network_regulators.tsv
│   ├── network_targets.tsv
│   ├── network_topology.tsv
│   └── network_source_counts.tsv
└── [visualization heatmaps & dotplots]
```

### Publication Figures
```
results/
├── volcano_combined.png/svg      # Volcano plots (RNA + Ribo)
├── volcano_rna.png/svg           # RNA-seq volcano
├── volcano_ribo.png/svg          # Ribo-seq volcano
├── scatter_translatome.png/svg   # RNA vs Ribo logFC scatter
├── curated_theme_*.png/svg       # Curated gene heatmaps
├── gene_logfc_heatmap.png/svg    # Gene-level logFC heatmap
├── *_heatmap*.png/svg            # Enrichment heatmaps
├── ridgeplot_*.png/svg           # Pathway ridge plots
├── summary_*.png                 # Publication summary figures
└── figure_captions.txt           # Manuscript figure captions
```

## Enrichment Analysis Details

### Over-Representation Analysis (ORA)
- **Method**: Hypergeometric test via clusterProfiler::enricher()
- **Test statistic**: P-value from hypergeometric distribution
- **Multiple testing correction**: Benjamini-Hochberg (FDR)
- **Significance threshold**: p.adjust < 0.1, gene set size ≥ 3

### Gene Sets
1. **GO terms**: Biological Process, Cellular Component, Molecular Function
2. **Pathways**: KEGG pathway definitions (Mycobrowser)
3. **Regulons**: Curated M. tuberculosis regulons (36 total)
4. **iModulons**: Systems biology modules from M. tuberculosis modulome (80 total)

### Gene Set Background
The background universe for ORA is the set of all genes with:
- Non-NA logFC and p-adjust values in the DE analysis
- Unique representation in annotation files

## Technical Notes

### Differential Expression
- **Tool**: DESeq2
- **Design**: Two-factor design (Condition: Control vs. Capreomycin)
- **Normalization**: DESeq2 median-of-ratios
- **Dispersion estimation**: Local fit
- **No manual filtering**: DESeq2's independent filtering automatically handles unreliable low-count genes

### Gene Categorization
After DE analysis, genes are classified as:
- **Co-regulated**: Significant in both RNA-seq AND Ribo-seq, same direction
- **Only RNA-seq**: Significant RNA-seq, not Ribo-seq
- **Only Ribo-seq**: Significant Ribo-seq, not RNA-seq
- **Opposite**: Significant in both but opposite directions
- **Not significant**: Doesn't meet significance threshold

### Visualization Standards
- **Color scheme**: 
  - Upregulated: Red (#FF6961)
  - Downregulated: Blue (#4285F4)
  - Not significant: Gray
- **Font sizes**: 24-30pt for publication quality
- **Resolution**: 300 dpi PNG, SVG for vector graphics

## Customization

### Change Sample Names
Edit `00_config.R`:
```r
rna_samples <- c("my_rna_1", "my_rna_2", ...)
ribo_samples <- c("my_ribo_1", "my_ribo_2", ...)
```

### Use Different Thresholds
```r
fdr_cutoff <- 0.05        # Stricter FDR
logfc_cutoff <- 1.0       # Higher fold-change requirement
```

### Add Custom Gene Sets
1. Create `.gmt` file in `pathways/`
2. Reference in enrichment script, e.g.:
```r
custom_geneset <- read.delim("pathways/my_custom_set.gmt")
```

### Run Subset of Steps
Comment out unwanted steps in `run_pipeline.R`:
```r
# run_step("04_upset_plot.R")  # Skip this step
# run_step("05b_volcano_plots.R")
run_step("05_enrichment_go_pathway.R")  # Run this
```

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| "File not found" | Wrong path or missing input | Check `00_config.R` paths and `data/` structure |
| "No significant genes" | Threshold too strict | Lower fdr_cutoff or logfc_cutoff |
| "No enrichment results" | Gene set too small after filtering | Check gene_map.tsv coverage, lower ora_min_gs_size |
| "Package not found" | Missing R package | See INSTALLATION.md for setup |
| "Out of memory" | Insufficient RAM | Use 16+ GB system or close other apps |
| "NaN in calculations" | Annotation file format issue | Verify .gmt, .csv, .tsv formatting |

## Pipeline Dependencies

**R packages** (version 4.0+):
- DESeq2 — Differential expression
- ggplot2, ggrepel, gridExtra, ggridges — Visualization
- dplyr, tidyr — Data manipulation
- clusterProfiler — Enrichment analysis
- ComplexHeatmap — Advanced heatmaps
- openxlsx — Excel export
- viridis — Color scales

**Data files** (in pathways/):
- gene_map.tsv — Gene ID mapping
- go_terms_genes.txt — GO annotations
- pathway_genes.txt — KEGG pathway annotations
- selected_regulon.gmt — Regulon definitions
- transcrnetwork_modulome.csv — Regulome network
- imodulons.tsv — iModulon mappings

## Related Documentation

- **../README.md** — Project overview and methods
- **../INSTALLATION.md** — Setup and installation
- **../QUICKSTART.md** — Quick reference guide
- **04_enrichment_common.R** — Shared ORA utilities (not a pipeline step)

---

**Pipeline version**: 1.0  
**Last updated**: 2026-10-28  
**Maintainer**: [Your name]
