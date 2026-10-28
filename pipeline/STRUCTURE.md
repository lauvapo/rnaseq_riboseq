# Pipeline Directory Structure

A clean, organized layout for the M. tuberculosis capreomycin analysis pipeline.

## Overview

```
pipeline/
├── 00_config.R              ⚙️  Configuration (shared by all steps)
├── run_pipeline.R           ▶️  Master driver (runs all steps)
├── STRUCTURE.md             📋 This file
├── README.md                📖 Technical documentation
│
├── utils/                   🔧 Shared utilities
│   └── 04_enrichment_common.R  (ORA helpers used by enrichment steps)
│
├── quality_control/         ✅ QC & validation
│   └── 02_quality_control_correlations.R
│
├── differential_expression/ 📊 DE analysis
│   └── 03_differential_expression.R
│
├── enrichment_analysis/     🧬 Gene set enrichment
│   ├── 05_transcription_network_analysis.R    (network topology)
│   ├── 05_enrichment_go_pathway.R             (GO/KEGG pathways)
│   ├── 05_enrichment_regulon.R                (regulons)
│   └── 05_enrichment_imodulons.R              (iModulons)
│
├── visualization/           📈 All plots & heatmaps
│   ├── 04_upset_plot.R                        (DE overlap)
│   ├── 05b_volcano_plots.R                    (RNA-seq vs Ribo-seq)
│   ├── 05c_scatter_transcriptome_translatome.R
│   ├── 05d_gene_category_heatmaps.R           (curated themes)
│   ├── 08_enrichment_heatmaps.R               (term × geneset)
│   ├── 09_enrichment_logfc_heatmaps.R         (mean logFC per term)
│   ├── 10_gene_logfc_heatmaps.R               (gene-level logFC)
│   ├── 11_curated_gene_heatmap.R              (curated genes)
│   └── 12_pathway_ridgeplots.R                (pathway distributions)
│
└── export_summary/          💾 Results export
    ├── 06_export_spreadsheets.R               (Excel workbook)
    └── 07_summary_figures.R                   (publication figures)
```

## Pipeline Execution Order

The master script `run_pipeline.R` runs these steps in order:

### 1. Quality Control (5 min)
```
quality_control/02_quality_control_correlations.R
```
- Replicate correlations
- RNA-seq vs Ribo-seq correlation
- Read length distributions

### 2. Differential Expression (10 min)
```
differential_expression/03_differential_expression.R
```
- DESeq2 analysis (RNA-seq)
- DESeq2 analysis (Ribo-seq)
- Gene classification (co-regulated, RNA-only, Ribo-only, opposite)

### 3. Visualization - Data Exploration (15 min)
```
visualization/04_upset_plot.R
visualization/05b_volcano_plots.R
visualization/05c_scatter_transcriptome_translatome.R
visualization/05d_gene_category_heatmaps.R
```
- UpSet plot of DE overlaps
- Volcano plots
- Transcriptome vs translatome scatter plot
- Curated gene heatmaps (5 biological themes)

### 4. Enrichment Analysis (15 min)
```
enrichment_analysis/05_transcription_network_analysis.R
enrichment_analysis/05_enrichment_go_pathway.R
enrichment_analysis/05_enrichment_regulon.R
enrichment_analysis/05_enrichment_imodulons.R
```
- Network topology analysis (regulome)
- Gene Ontology enrichment (BP/CC/MF)
- KEGG pathway enrichment
- Regulon-based enrichment (36 regulons)
- iModulon enrichment (80 iModulons)

### 5. Visualization - Enrichment Results (10 min)
```
visualization/08_enrichment_heatmaps.R
visualization/09_enrichment_logfc_heatmaps.R
visualization/10_gene_logfc_heatmaps.R
visualization/11_curated_gene_heatmap.R
visualization/12_pathway_ridgeplots.R
```
- Term × gene-set ORA heatmap
- Mean logFC heatmaps per term
- Gene-level logFC heatmap
- Curated gene heatmap
- Pathway ridge plots

### 6. Export & Summary (5 min)
```
export_summary/06_export_spreadsheets.R
export_summary/07_summary_figures.R
```
- Excel workbook with all results
- Publication-ready summary figures

**Total runtime**: 30-60 minutes (depending on system)

## File Organization Rationale

### Why organize this way?

1. **Logical grouping**: Scripts grouped by analysis purpose, not execution order
2. **Easier navigation**: Find visualization scripts in one place, enrichment in another
3. **Scalability**: Easy to add new scripts without cluttering the root folder
4. **Clarity**: Directory names describe what happens inside
5. **Maintainability**: Simpler to review and update specific analysis modules

### Directory Purposes

| Folder | Purpose | # Scripts |
|--------|---------|-----------|
| utils/ | Shared functions, helpers, constants | 1 |
| quality_control/ | Data validation, QC metrics | 1 |
| differential_expression/ | Statistical DE analysis (DESeq2) | 1 |
| enrichment_analysis/ | Gene set enrichment (ORA) | 4 |
| visualization/ | All plots, heatmaps, figures | 9 |
| export_summary/ | Results export, manuscript figures | 2 |

## Updating the Pipeline

### Add a new visualization
1. Create `visualization/05x_new_plot.R`
2. Add to `run_pipeline.R`: `run_step("visualization/05x_new_plot.R")`
3. Update this file with description

### Add a new enrichment type
1. Create `enrichment_analysis/05_enrichment_newtype.R`
2. Source utils: `source(file.path(pipeline_dir, "utils/04_enrichment_common.R"))`
3. Add to `run_pipeline.R` and export script

### Share code across scripts
1. Add utility functions to `utils/04_enrichment_common.R`
2. Source from other scripts: `source(file.path(pipeline_dir, "utils/04_enrichment_common.R"))`

## Path Resolution

All scripts use `pipeline_dir` variable (set in `00_config.R`) to find other scripts and resources:

```r
# In any subscript:
source(file.path(pipeline_dir, "utils/04_enrichment_common.R"))
source(file.path(pipeline_dir, "quality_control/02_quality_control_correlations.R"))

# This works regardless of where run_pipeline.R is invoked from
```

## Deleted Files

These files were removed as they are no longer used:

- **05_enrichment_network.R** — Network enrichment (removed from pipeline)
- **05e_multipanel_heatmaps.R** — Multipanel heatmaps (not called by pipeline)

Legacy code is in the `remove/deprecated/` folder if needed for reference.

---

**Version**: 1.0  
**Last updated**: 2026-10-28
