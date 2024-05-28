#!/bin/bash

# Diretórios e variáveis
FASTQ_DIR="/path/to/fastq_files"
FASTQC_DIR="/path/to/Fastqc_analysis"
TRIMMOMATIC_DIR="/path/to/trimmomatic"
HISAT2_INDEX_DIR="/path/to/Index_Hisat2"
REFERENCE_GENOME="/path/to/genome_reference"
GFF3_FILE="/path/to/Gff3_file"
OUTPUT_DIR="/path/to//output"
GTF_FILE="/path/to/annotation.gtf"
SCRIPTS_DIR="/path/to/scripts"

# Verificar se os diretórios de saída existem, caso contrário, criá-los
mkdir -p $FASTQC_DIR
mkdir -p $OUTPUT_DIR/trimmed
mkdir -p $OUTPUT_DIR/hisat2
mkdir -p $OUTPUT_DIR/samtools
mkdir -p $OUTPUT_DIR/prepDE

# FastQC
echo "Running FastQC..."
fastqc -t 15 -o $FASTQC_DIR $FASTQ_DIR/*.fastq.gz

# MultiQC
echo "Running MultiQC..."
multiqc $FASTQC_DIR -o $FASTQC_DIR

# Trimmomatic
echo "Running Trimmomatic..."
for sample in $FASTQ_DIR/*.fastq.gz; do
    base=$(basename $sample ".fastq.gz")
    java -jar $TRIMMOMATIC_DIR/trimmomatic-0.36.jar PE -threads 15 -phred64 \
    $FASTQ_DIR/${base}_R1_001.fastq.gz $FASTQ_DIR/${base}_R2_001.fastq.gz \
    $OUTPUT_DIR/trimmed/${base}_R1_paired.fastq.gz $OUTPUT_DIR/trimmed/${base}_R1_unpaired.fastq.gz \
    $OUTPUT_DIR/trimmed/${base}_R2_paired.fastq.gz $OUTPUT_DIR/trimmed/${base}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 MINLEN:36
done

# Indexação do genoma de referência com HISAT2
echo "Indexing reference genome with HISAT2..."
hisat2_extract_splice_sites.py $HISAT2_INDEX_DIR/$GFF3_FILE > $HISAT2_INDEX_DIR/splicesites.tsv
hisat2_extract_exons.py $HISAT2_INDEX_DIR/$GFF3_FILE > $HISAT2_INDEX_DIR/exons.tsv
hisat2-build -p 20 --ss $HISAT2_INDEX_DIR/splicesites.tsv --exon $HISAT2_INDEX_DIR/exons.tsv $HISAT2_INDEX_DIR/$REFERENCE_GENOME $HISAT2_INDEX_DIR/index

# HISAT2
echo "Running HISAT2..."
for sample in $OUTPUT_DIR/trimmed/*_paired.fastq.gz; do
    base=$(basename $sample "_paired.fastq.gz")
    hisat2 -p 15 -x $HISAT2_INDEX_DIR/index \
    -1 $OUTPUT_DIR/trimmed/${base}_R1_paired.fastq.gz \
    -2 $OUTPUT_DIR/trimmed/${base}_R2_paired.fastq.gz \
    -S $OUTPUT_DIR/hisat2/${base}.sam
done

# Samtools
echo "Running Samtools..."
for samfile in $OUTPUT_DIR/hisat2/*.sam; do
    base=$(basename $samfile ".sam")
    samtools view -bS $samfile > $OUTPUT_DIR/samtools/${base}.bam
    samtools sort $OUTPUT_DIR/samtools/${base}.bam -o $OUTPUT_DIR/samtools/${base}_sorted.bam
    samtools index $OUTPUT_DIR/samtools/${base}_sorted.bam
done

# StringTie
echo "Running StringTie..."
for bamfile in $OUTPUT_DIR/samtools/*_sorted.bam; do
    base=$(basename $bamfile "_sorted.bam")
    stringtie -e -B -p 60 -G $HISAT2_INDEX_DIR/$GFF3_FILE \
    -o $OUTPUT_DIR/prepDE/${base}.gtf $bamfile
done

# prepDE.py
echo "Running prepDE.py..."
python $SCRIPTS_DIR/prepDE.py -i $OUTPUT_DIR/samtools -g $OUTPUT_DIR/prepDE/gene_count_matrix.csv -t $OUTPUT_DIR/prepDE/transcript_count_matrix.csv

echo "Pipeline finished!"
