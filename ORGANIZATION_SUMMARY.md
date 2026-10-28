# Project Organization Summary

This document describes the cleanup and reorganization performed on the capreomycin remapping project.

## Changes Made

### 1. **Consolidated Input Data**
   - Moved: `RNA_seq/` and `Ribo_seq/` → `data/`
   - Purpose: Separate raw inputs from analysis outputs
   - Updated: Pipeline paths in `pipeline/00_config.R`

### 2. **Organized Pipeline Outputs**
   - Moved: `significative_genes/` and `enrichment_results/` → `results/`
   - Purpose: All pipeline results now live in a single output directory
   - Updated: Pipeline config to write to `results/significative_genes/` and `results/enrichment_results/`

### 3. **Grouped Exploratory Analyses**
   - Moved: All `.Rmd` notebooks → `analyses/`
   - Moved: `plots/`, `heatmaps/`, `images/`, `pagerank/` (notebook outputs) → `analyses/outputs/`
   - Purpose: Separate exploratory work from main pipeline

### 4. **Organized Utilities**
   - Created: `scripts/` directory
   - Moved: `homoplasy_counter.py` → `scripts/`
   - Purpose: Central location for standalone analysis tools

### 5. **Organized Publication Materials**
   - Created: `supplementary/` directory
   - Moved: `final_figures/` → `supplementary/`
   - Purpose: Publication-ready outputs separate from exploratory work

### 6. **Cleaned Up Duplicates & Backups**
   - Moved: `enrichment_results (Copiar)`, `results (Copiar)`, `prueba/` → `backups/`
   - Deleted: `.Rhistory` and `.RData` files (session-specific clutter)
   - Purpose: Archive old work without cluttering the active project

### 7. **Consolidated Reference Data**
   - Existing: `pathways/` remains at root
   - Contains: GO terms, pathway annotations, regulatory network
   - These are inputs to the pipeline and remain unchanged

### 8. **Deprecated Analyses**
   - Existing: `deprecated/legacy_pagerank_outputs/` preserved
   - Purpose: Historical record of previous methodology changes

## Final Directory Structure

```
paper_capreomycin_remapping/
├── README.md                          # Project overview
├── .gitignore                         # Git exclusions
├── ORGANIZATION_SUMMARY.md            # This file
├── .claude/                           # Claude Code configuration
│
├── pipeline/                          # Main automated analysis
│   ├── 00_config.R                    # Parameters & paths (UPDATED)
│   ├── 01-11_*.R                      # Analysis scripts
│   ├── run_pipeline.R                 # Entry point
│   └── README.md                      # Pipeline documentation (UPDATED)
│
├── results/                           # PRIMARY OUTPUTS (pipeline)
│   ├── significative_genes/           # Gene significance tables
│   ├── enrichment_results/            # GO/pathway ORA results
│   ├── results/                       # Summary tables & plots
│   └── README.md                      # Results guide (NEW)
│
├── data/                              # Raw input data (REORGANIZED)
│   ├── RNA_seq/
│   │   ├── counts/
│   │   └── read_lengths/
│   ├── Ribo_seq/
│   │   ├── counts/
│   │   └── read_lengths/
│   └── README.md                      # Data documentation (NEW)
│
├── analyses/                          # Exploratory notebooks (ORGANIZED)
│   ├── *.Rmd                          # Moved from root
│   ├── outputs/                       # Notebook outputs (ORGANIZED)
│   │   ├── plots/
│   │   ├── heatmaps/
│   │   ├── images/
│   │   └── pagerank/
│   └── README.md                      # Analyses guide (NEW)
│
├── supplementary/                     # Publication figures (NEW)
│   └── final_figures/
│       ├── 0.05/
│       └── 0.1/
│
├── scripts/                           # Utility scripts (NEW)
│   ├── homoplasy_counter.py           # Moved from root
│   └── README.md                      # Scripts guide (NEW)
│
├── pathways/                          # Reference data (unchanged)
│   ├── go_terms_genes.txt
│   ├── pathway_genes.txt
│   └── transcrnetwork_bustad.csv
│
├── deprecated/                        # Superseded analyses (unchanged)
│   └── legacy_pagerank_outputs/
│
└── backups/                           # Archived versions (NEW)
    ├── enrichment_results (Copiar)/
    ├── results (Copiar)/
    └── prueba/
```

## Migration Guide

### For existing R scripts
If you have custom scripts that reference the old paths, update:
- `RNA_seq/` → `data/RNA_seq/`
- `Ribo_seq/` → `data/Ribo_seq/`
- `significative_genes/` → `results/significative_genes/`
- `enrichment_results/` → `results/enrichment_results/`

The pipeline itself (`pipeline/00_config.R`) has already been updated.

### For running the pipeline
No changes needed to the command:
```bash
Rscript pipeline/run_pipeline.R
```
All paths are configured in `pipeline/00_config.R`.

### For exploratory notebooks
Notebooks in `analyses/` should still work as-is (they load paths from `../pipeline/00_config.R`). If a notebook was previously writing to `plots/`, it will now write to `analyses/outputs/plots/` (adjust notebook `output_dir` settings if needed).

## Benefits of This Organization

1. **Clear separation of concerns**: Inputs, pipeline, outputs, explorations, utilities
2. **Easier navigation**: Everyone knows where to find results, data, and reference files
3. **Cleaner backups**: Old versions archived in one place, not scattered at root
4. **Reproducible paths**: Single source of truth for file paths (`pipeline/00_config.R`)
5. **Git-friendly**: `.gitignore` configured to exclude session files and large outputs
6. **Scalable**: Easy to add new notebooks, scripts, or analyses without cluttering root

## Next Steps

1. Verify the pipeline runs correctly with the new paths:
   ```bash
   cd pipeline && Rscript run_pipeline.R
   ```

2. Test that exploratory notebooks can still access their data (they should, since they source `../pipeline/00_config.R`)

3. Update any external scripts or documentation that reference the old file structure

4. (Optional) Commit this reorganization to git:
   ```bash
   git add -A && git commit -m "refactor: reorganize project structure for clarity"
   ```
