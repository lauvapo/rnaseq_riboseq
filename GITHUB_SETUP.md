# GitHub Setup Checklist

Your pipeline is now GitHub-ready! Follow these steps to publish it.

## Before Publishing to GitHub

- [ ] Review and update `README.md` with your specific details:
  - [ ] Your name and affiliation
  - [ ] Citation information (once paper is published)
  - [ ] Contact email

- [ ] Update `INSTALLATION.md`:
  - [ ] Contact information for support

- [ ] Update `LICENSE` with your name and year

- [ ] Verify `.gitignore` covers all sensitive/large files:
  - [ ] `data/` (raw counts) ✓
  - [ ] `results/` (outputs) ✓
  - [ ] `.Rhistory`, `Rplots.pdf` ✓

## Initialize Git Repository

```bash
cd /path/to/paper_capreomycin_remapping

# Initialize git
git init

# Add all tracking files (respects .gitignore)
git add .

# Check what will be committed
git status

# Create initial commit
git commit -m "Initial commit: M. tuberculosis capreomycin RNA-seq/Ribo-seq pipeline"
```

## Create GitHub Repository

1. **Go to GitHub** → [github.com/new](https://github.com/new)
2. **Repository name**: `capreomycin-pipeline` (or your preference)
3. **Description**: "RNA-seq and Ribo-seq analysis pipeline for capreomycin-induced transcriptomic changes in M. tuberculosis H37Rv"
4. **Visibility**: Public
5. **Do NOT initialize** with README (we already have one)
6. Click **Create repository**

## Push to GitHub

```bash
# Add remote
git remote add origin https://github.com/YOUR_USERNAME/capreomycin-pipeline.git

# Push
git branch -M main
git push -u origin main
```

## Add Topics (GitHub)

In your repository settings, add these topics:
- `bioinformatics`
- `transcriptomics`
- `rna-seq`
- `ribosome-profiling`
- `mycobacterium`
- `r`
- `reproducible-research`

## Enable GitHub Pages (Optional)

To host pipeline documentation:
1. Go to **Settings** → **Pages**
2. Select **Branch**: `main`
3. Select **Folder**: `/docs` (optional)
4. Click **Save**

Your documentation will be at: `https://YOUR_USERNAME.github.io/capreomycin-pipeline/`

## Create Release (Optional)

```bash
# Tag the release
git tag -a v1.0 -m "First stable release"

# Push tags
git push origin v1.0

# On GitHub, go to Releases → Draft a new release → select v1.0
# Add release notes and changelog
```

## File Structure for GitHub

**What's Included** (uploaded to GitHub):
```
capreomycin-pipeline/
├── .gitignore
├── LICENSE
├── README.md
├── INSTALLATION.md
├── QUICKSTART.md
├── GITHUB_SETUP.md (this file)
├── pipeline/
│   ├── 00_config.R
│   ├── run_pipeline.R
│   └── [all other analysis scripts]
├── pathways/
│   ├── gene_map.tsv
│   ├── go_terms_genes.txt
│   ├── pathway_genes.txt
│   ├── selected_regulon.gmt
│   ├── transcrnetwork_modulome.csv
│   └── imodulons.tsv
└── scripts/
    └── [utility scripts]
```

**What's NOT Included** (.gitignore):
```
data/                      # Raw RNA-seq/Ribo-seq counts
results/                   # Pipeline outputs
remove/                    # Legacy/deprecated files
.Rhistory, Rplots.pdf     # R artifacts
.claude/settings.local.json # Local settings
```

## Sharing with Collaborators

### Grant Access (Private Repo)
1. **Settings** → **Collaborators** → **Add people**
2. Choose permission level (Write for active contributors)

### Using Branches for Development
```bash
# Create feature branch
git checkout -b feature/new-analysis

# Make changes, commit
git add .
git commit -m "Add new enrichment method"

# Push branch
git push origin feature/new-analysis

# Create Pull Request on GitHub
```

## Continuous Integration (Optional)

To automatically check code quality on every push:

Create `.github/workflows/check.yml`:
```yaml
name: Check Pipeline

on: [push, pull_request]

jobs:
  check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: r-lib/actions/setup-r@v2
      - name: Install dependencies
        run: Rscript -e "install.packages(c('ggplot2', 'DESeq2'))"
      - name: Syntax check
        run: Rscript -e "library(lintr); lint_dir('pipeline')"
```

## Sharing Results (Data Only)

If you want to share results without exposing raw data:

1. **Zenodo** (recommended for datasets):
   - Upload `results/DE_and_enrichment_results.xlsx`
   - Get DOI for citation

2. **Open Science Framework (OSF)**:
   - Host raw data + results
   - Link from GitHub README

3. **GitHub Releases**:
   - Attach Excel workbook to release

## Documentation Best Practices

✓ Keep README updated with any changes  
✓ Add section to README for each new feature  
✓ Document custom parameters in INSTALLATION.md  
✓ Include example outputs  
✓ Link to published paper once available  

## SEO & Discoverability

To improve findability:

1. **GitHub Topics** (added above)
2. **Keywords in README** (bioinformatics, RNA-seq, etc.)
3. **Link from your institution's website**
4. **Submit to bioinformatics portals**:
   - [Bioconda](https://bioconda.github.io/) (optional)
   - [Omictools.com](http://omictools.com/)
   - [Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/)

## Maintenance

### Regular Updates
- [ ] Keep R packages updated
- [ ] Test pipeline quarterly
- [ ] Update docs with new features
- [ ] Fix reported issues promptly

### Version Numbering
Follow [Semantic Versioning](https://semver.org/):
- **v1.0.0** — First stable release
- **v1.1.0** — New features added
- **v1.0.1** — Bug fixes only
- **v2.0.0** — Major breaking changes

## Example GitHub Profile

```markdown
## Capreomycin Pipeline

**Reproducible RNA-seq & Ribo-seq analysis for M. tuberculosis**

- [📖 Documentation](https://github.com/yourusername/capreomycin-pipeline)
- [🚀 Quick Start](https://github.com/yourusername/capreomycin-pipeline/blob/main/QUICKSTART.md)
- [📥 Results](https://zenodo.org/record/xxxxx)

### Citation
[Your paper citation once published]
```

## Questions?

- **Setup issues**: Check GitHub's [Getting Started](https://docs.github.com/en/get-started)
- **R package issues**: See INSTALLATION.md
- **Data format issues**: See QUICKSTART.md

---

**Your pipeline is ready to share with the world!** 🚀

Once published, you can link to it from your CV, institution website, and paper supplementary materials.
