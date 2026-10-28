#!/bin/bash

PROCS=3
DIR=.
SAMPLES=$(ls ./fastq/*_R1_001.fastq | xargs -n 1 basename | sed 's/_R1_001.fastq//' | sort -u)
REF=${DIR}/genome_resources/Mycobacterium_tuberculosis_H37Rv_
SAMTOOLS=/data/ThePipeline_v2/Programs/samtools-1.15/samtools

echo ">>>>>>>>>>>>> STARTING PREPROCESSING <<<<<<<<<<<<<<<"
mkdir -p ${DIR}
mkdir -p ${DIR}/alignments
mkdir -p ${DIR}/counts
mkdir -p ${DIR}/genome_resources
mkdir -p ${DIR}/fastq
mkdir -p ${DIR}/read_lengths

cd ${DIR}/genome_resources

if [[ ! -f Mycobacterium_tuberculosis_H37Rv_genome_v4.fasta ]]; then
    echo ">>>>>>>>>>>>> Downloading Mycobacterium tuberculosis H37Rv genome v4 <<<<<<<<<<<<<<<"
    wget -O Mycobacterium_tuberculosis_H37Rv_genome_v4.fasta \
    'https://mycobrowser.epfl.ch/releases/4/get_file?dir=fasta&file=Mycobacterium_tuberculosis_H37Rv.fasta'
fi

echo ">>>>>>>>>>>>> Building indexes for Mycobacterium_tuberculosis_H37Rv_genome_v4 genome <<<<<<<<<<<<<<<"
bowtie2-build Mycobacterium_tuberculosis_H37Rv_genome_v4.fasta Mycobacterium_tuberculosis_H37Rv_genome_v4.fasta
bowtie2-build Mycobacterium_tuberculosis_H37Rv_tRNA_rRNA.fasta Mycobacterium_tuberculosis_H37Rv_tRNA_rRNA.fasta
echo ">>>>>>>>>>>>> Indexes built <<<<<<<<<<<<<<<"

cd ..

    
for SAMPLE in $SAMPLES
do
    echo ">>>>>>>>>>>>> Analyzing sample '$SAMPLE' <<<<<<<<<<<<<<<"
    echo ">>>>>> Trimming adapters with cutadapt..."

    INFASTQ=./fastq/${SAMPLE}_R1_001.fastq

    nice -n 5 cutadapt \
        -a AAGATCGGAAGAGCACACGTCT \
        -O 6 \
        -m 10 \
        -e 0.15 \
        -o ${DIR}/fastq/${SAMPLE}_R1_001_trimmed.fastq.gz \
        --too-short-output ${DIR}/fastq/${SAMPLE}_R1_001_length_discarded.fastq.gz \
        --untrimmed-output ${DIR}/fastq/${SAMPLE}_R1_001_untrimmed.fastq.gz \
        "$INFASTQ"

    echo ">>>>>> Aligning tRNA-rRNA reads with bowtie2..."
    nice -n 5 bowtie2\
    		-x ${REF}tRNA_rRNA.fasta\
    		--sensitive\
    		--al-gz ${DIR}/alignments/${SAMPLE}.tRNArRNA.fastq.gz\
        	--un-gz ${DIR}/alignments/${SAMPLE}.ncRNAmRNA.fastq.gz\
        	-U ${DIR}/fastq/${SAMPLE}_R1_001_trimmed.fastq.gz\
        	-q\
        	-p $PROCS\
        	-S ${DIR}/alignments/${SAMPLE}.tRNArRNA.sam
    
    echo ">>>>>> Converting tRNA-rRNA SAM to BAM with samtools..."
    nice -n 5 $SAMTOOLS view\
    		-b\
    		-S\
    		-@ $PROCS\
    		-o ${DIR}/alignments/${SAMPLE}.tRNArRNA.bam\
    		${DIR}/alignments/${SAMPLE}.tRNArRNA.sam
    
    rm ${DIR}/alignments/${SAMPLE}.tRNArRNA.sam

    echo ">>>>>> Aligning ncRNA-mRNA reads with bowtie2..."
    nice -n 5 bowtie2\
    		-x ${REF}genome_v4.fasta\
    		-U ${DIR}/alignments/${SAMPLE}.ncRNAmRNA.fastq.gz\
        	-N 1\
        	-q\
        	-p $PROCS\
        	-S ${DIR}/alignments/${SAMPLE}.sam
        	
    echo ">>>>>> Converting ncRNA-mRNA SAM to BAM and filtering by length with samtools..."
    nice -n 5 $SAMTOOLS view\
    		-e "qlen>11 && mapq>10"\
    		-b\
    		-S\
    		-@ $PROCS\
    		-o ${DIR}/alignments/${SAMPLE}.bam\
    		${DIR}/alignments/${SAMPLE}.sam
    
    rm ${DIR}/alignments/${SAMPLE}.sam

    echo ">>>>>> Sorting BAM with samtools"
    nice -n 5 $SAMTOOLS sort\
    		-n\
    		-@ $PROCS\
    		-o ${DIR}/alignments/${SAMPLE}.sortname.bam\
    		${DIR}/alignments/${SAMPLE}.bam
    
    rm ${DIR}/alignments/${SAMPLE}.bam
    nice -n 5 featureCounts -T $PROCS -t CDS -g Name -s 1 --minOverlap 1 -O -M --fraction -a ${REF}v3.gff \
        -o ${DIR}/counts/${SAMPLE}.dedup.csv ${DIR}/alignments/${SAMPLE}.sortname.bam
    nice -n 5 featureCounts -T $PROCS -t five_prime_UTR -g Name -s 1 --minOverlap 1 -O -M --fraction -a ${DIR}/genome_resources/5UTR_annotation_2.7.gff \
        -o ${DIR}/counts/${SAMPLE}_5UTR.dedup.csv ${DIR}/alignments/${SAMPLE}.sortname.bam
    nice -n 5 featureCounts -T $PROCS -t ncRNA -g Name -s 1 --minOverlap 1 -O -M --fraction -a ${DIR}/genome_resources/ncRNA_annotation_2.7.gff \
        -o ${DIR}/counts/${SAMPLE}_ncRNA.dedup.csv ${DIR}/alignments/${SAMPLE}.sortname.bam
    nice -n 5 featureCounts -T $PROCS -t intergenic -g Name -s 1 --minOverlap 1 -O -M --fraction -a ${DIR}/genome_resources/IGR_annotation_2.7_stranded.gff \
        -o ${DIR}/counts/${SAMPLE}_IGR.dedup.csv ${DIR}/alignments/${SAMPLE}.sortname.bam
    nice -n 5 featureCounts -T $PROCS -t CDS -g Name -s 1 --minOverlap 1 -O -M --fraction -a ${DIR}/genome_resources/CDS_annotation_2.7.gff \
        -o ${DIR}/counts/${SAMPLE}_CDS.dedup.csv ${DIR}/alignments/${SAMPLE}.sortname.bam
    echo ">>>>>>>>>>>>> Sample '$SAMPLE' analyzed <<<<<<<<<<<<<<<"	
    
done


for bam in ${DIR}/alignments/*sortname.bam; do
    sample=$(basename "$bam" .bam)
    
    echo ">>>>>>>>>>>>> Calculating read lengths from $sample <<<<<<<<<<<<<<<"

    # 1) Length (format: length  count)
    samtools view -F 4 "$bam" | \
    awk '{print length($10)}' | \
    sort -n | uniq -c | \
    awk '{print $2, $1}' > read_lengths/${sample}.lengths.txt

    # 2) Length (format: count  length)
    samtools view -F 4 "$bam" | \
    awk '{print length($10)}' | \
    sort -n | uniq -c > read_lengths/${sample}_mappedread-length-distribution.txt
    
done

echo ">>>>>>>>>>>>> PREPROCESSING FINISHED <<<<<<<<<<<<<<<"

