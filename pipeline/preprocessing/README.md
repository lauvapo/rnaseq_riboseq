# Raw Data Preprocessing Scripts

Scripts for processing raw FASTQ reads into count tables for differential expression analysis.

## Overview

These scripts handle the complete preprocessing pipeline:

1. **Adapter trimming** (cutadapt)
2. **Quality control** (FastQC)
3. **rRNA/tRNA filtering** (bowtie2)
4. **Genome alignment** (bowtie2)
5. **SAM to BAM conversion** (samtools)
6. **PCR duplicate removal** (samtools)
7. **Feature counting** (featureCounts)

## Prerequisites

### Conda Environment (Recommended)

```bash
# From project root
conda env create -f environment.yml
conda activate capreomycin-pipeline
```

This installs:
- bowtie2 (alignment)
- samtools (BAM manipulation)
- cutadapt (adapter trimming)
- subread/featureCounts (counting)
- fastqc (quality control)

### Manual Installation (if not using Conda)

```bash
# Install each tool individually
conda install -c bioconda bowtie2
conda install -c bioconda samtools
conda install -c bioconda cutadapt
conda install -c bioconda subread
conda install -c bioconda fastqc
```

## Input Requirements

### Directory Structure

```
project_root/
├── fastq/
│   ├── sample1_R1_001.fastq.gz
│   ├── sample1_R2_001.fastq.gz (if paired-end)
│   ├── sample2_R1_001.fastq.gz
│   └── sample2_R2_001.fastq.gz
└── genome_resources/
    ├── Mycobacterium_tuberculosis_H37Rv_genome_v4.fasta
    └── Mycobacterium_tuberculosis_H37Rv_tRNA_rRNA.fasta
```

### FASTQ File Naming

Standard Illumina format:
- `SAMPLENAME_R1_001.fastq.gz` (read 1, forward)
- `SAMPLENAME_R2_001.fastq.gz` (read 2, reverse, optional)

The scripts automatically detect samples from this naming convention.

## Running Preprocessing

### 1. RNA-seq Preprocessing

```bash
# Activate environment
conda activate capreomycin-pipeline

# Run from project root (where fastq/ directory exists)
bash pipeline/preprocessing/process_fastq_rnaseq.sh
```

**Output:**
- `fastq/*_trimmed.fastq.gz` — Trimmed reads
- `alignments/*_genome.bam` — Genome-aligned BAM files
- `counts/*_CDS.dedup.csv` — Feature counts (input for DESeq2)

### 2. Ribo-seq Preprocessing

```bash
bash pipeline/preprocessing/process_fastq_riboseq.sh
```

Same output structure, generated in separate directories.

## Script Details

### process_fastq_rnaseq.sh

Processes total RNA-seq reads:

1. **Adapter trimming** with cutadapt
   - Removes Illumina adapter (AAGATCGGAAGAGCACACGTCT)
   - Minimum read length: 10 bp
   - Minimum overlap: 6 bp

2. **rRNA/tRNA filtering**
   - Aligns to rRNA/tRNA sequences with bowtie2
   - Removes contaminating ribosomal reads

3. **Genome alignment**
   - Aligns to M. tuberculosis H37Rv reference genome
   - Local alignment mode (--local)
   - Removes duplicates with samtools

4. **Feature counting**
   - Uses featureCounts from subread
   - Counts reads per CDS (protein-coding gene)
   - Generates count matrix for DE analysis

### process_fastq_riboseq.sh

Processes ribosome profiling (Ribo-seq) reads:

- Similar workflow to RNA-seq
- Parameters optimized for footprint analysis
- Generates separate count tables

## Output Files

After preprocessing, you'll have:

```
alignments/
├── SAMPLE1_genome.bam     ← Deduplicated, genome-aligned reads
├── SAMPLE1_genome.bam.bai
├── SAMPLE2_genome.bam
└── ...

counts/
├── SAMPLE1_CDS.dedup.csv  ← Feature counts (→ DESeq2)
├── SAMPLE2_CDS.dedup.csv
└── ...

read_lengths/
├── SAMPLE1_lengths.txt    ← For QC analysis
└── ...
```

**The `*_CDS.dedup.csv` files** are what the main pipeline expects.

## Integration with Main Pipeline

### 1. Preprocess raw FASTQ files
```bash
bash pipeline/preprocessing/process_fastq_rnaseq.sh
bash pipeline/preprocessing/process_fastq_riboseq.sh
```

### 2. Copy counts to expected locations
```bash
# RNA-seq counts
cp counts/RNA/*_CDS.dedup.csv data/RNA_seq/counts/CDS/

# Ribo-seq counts
cp counts/Ribo/*_CDS.dedup.csv data/Ribo_seq/counts/CDS/
```

### 3. Run main analysis pipeline
```bash
Rscript pipeline/run_pipeline.R
```

## Customization

### Change trimming parameters

Edit the cutadapt command in the script:

```bash
cutadapt \
    -a ADAPTER_SEQUENCE \    # Your adapter (default: Illumina)
    -O 6 \                   # Min overlap length
    -m 10 \                  # Min final read length
    -e 0.15 \                # Max error rate
    -j 0 \                   # Threads (0 = all)
    -o output.fastq.gz \
    input.fastq.gz
```

### Change alignment parameters

Edit the bowtie2 command:

```bash
bowtie2 \
    -x GENOME_INDEX \
    -U reads.fastq.gz \
    --local \              # Local vs end-to-end alignment
    -X 1000 \             # Max fragment length
    -p 4                  # Threads
```

### Change counting parameters

Edit the featureCounts command:

```bash
featureCounts \
    -a ANNOTATION.gff \
    -o output.csv \
    -t exon \            # Feature type
    -g gene_id \         # Gene ID attribute
    -s 0 \               # Strandedness (0=unstranded)
    sample.bam
```

## Performance Optimization

### Parallel Processing

Both scripts use multiple threads. Adjust:

```bash
PROCS=8  # Default is 3, increase for faster processing
```

At the top of each script.

### Speed Benchmarks

On a 4-core machine with 8GB RAM:
- RNA-seq (40M reads): ~30-45 minutes
- Ribo-seq (30M reads): ~20-30 minutes

## Troubleshooting

### "bowtie2: command not found"
```bash
# Activate Conda environment
conda activate capreomycin-pipeline
```

### "No such file or directory: genome_resources/"
- Create the directory: `mkdir -p genome_resources`
- Script will auto-download genome reference

### "Permission denied" when running script
```bash
chmod +x pipeline/preprocessing/process_fastq_rnaseq.sh
chmod +x pipeline/preprocessing/process_fastq_riboseq.sh
```

### Out of memory
- Reduce `PROCS` (threads) at top of script
- Process samples one at a time

### Very slow processing
- Check CPU usage: `top` (should show ~100% × PROCS)
- Check disk I/O: `iostat` (if slow, move to SSD)
- Reduce `PROCS` if running other programs

## Reference Genome

Automatically downloaded from Mycobrowser:
- **Genome**: M. tuberculosis H37Rv v4 (FASTA)
- **Source**: https://mycobrowser.epfl.ch/
- **Size**: ~4.4 MB (uncompressed)

To use a different genome:

1. Edit script and replace download URL
2. Update bowtie2-build index name
3. Update featureCounts annotation file

## Output Format

The count tables have this format:

```
Geneid,Length,J_S4_CDS.dedup.bam
Rv0001,1309,145
Rv0002,1032,89
...
```

- **Column 1**: Gene locus tag (e.g., Rv0001)
- **Column 2**: Gene length (bp)
- **Column 3+**: Count per sample

This matches the format expected by the main pipeline (step 03_differential_expression.R).

## Quality Control

### Check alignment quality

```bash
samtools flagstat alignments/SAMPLE_genome.bam
```

Shows:
- % mapped reads
- % paired correctly (if applicable)
- % duplicates

### Check read length distribution

```bash
# Generated by script automatically
cat read_lengths/SAMPLE_lengths.txt | head -20
```

### Run FastQC on trimmed reads

```bash
conda activate capreomycin-pipeline
fastqc fastq/*_trimmed.fastq.gz
multiqc fastq/
```

Opens `multiqc_report.html` with comprehensive QC report.

## Citation

If using these scripts, cite the tools:

- Langmead & Salzberg (2012). Bowtie2. Nature Methods.
- Martin (2011). Cutadapt. EMBnet.journal.
- Li et al. (2009). Samtools. Bioinformatics.
- Liao et al. (2014). FeatureCounts. Bioinformatics.

---

**Next Step**: Run preprocessing, then proceed to main pipeline:
```bash
bash pipeline/preprocessing/process_fastq_rnaseq.sh
bash pipeline/preprocessing/process_fastq_riboseq.sh
Rscript pipeline/run_pipeline.R
```
