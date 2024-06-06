#!/bin/bash

# Directories and variables
INPUT_DIR="/path/to/RNAseq/files"
TRIMMOMATIC_DIR="${OUTPUT_DIR}/trimmomatic"
HISAT2_INDEX_DIR="/path/to/genome/reference"
REFERENCE_GENOME="genome.fna"
GFF3_FILE="annotation.gff"
OUTPUT_DIR="/path/to/results"
SCRIPTS_DIR="/path/to/scripts"

# Check if output directories exist, create them if they don't
mkdir -p $OUTPUT_DIR/fastqc
mkdir -p $OUTPUT_DIR/trimmed
mkdir -p $OUTPUT_DIR/hisat2
mkdir -p $OUTPUT_DIR/samtools
mkdir -p $OUTPUT_DIR/prepDE

# FastQC
#echo "Running FastQC..."
#fastqc -t 15 -o $OUTPUT_DIR/fastqc $INPUT_DIR/*.fastq

# Trimmomatic
echo "Running Trimmomatic..."
for sample in $INPUT_DIR/*_R1.fastq; do
    base=$(basename $sample "_R1.fastq")  # Get the base part of the filename
    
    # Define the corresponding file with suffix _R2.fastq
    sample_R2="${INPUT_DIR}/${base}_R2.fastq"
    
    java -jar /path/to/Trimmomatic-0.39/trimmomatic.jar PE -threads 15 \
    $sample $sample_R2 \
    $OUTPUT_DIR/trimmed/${base}_R1_paired.fastq $OUTPUT_DIR/trimmed/${base}_R1_unpaired.fastq \
    $OUTPUT_DIR/trimmed/${base}_R2_paired.fastq $OUTPUT_DIR/trimmed/${base}_R2_unpaired.fastq \
    ILLUMINACLIP:/path/to/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:20 TOPHRED33
done

# Index the reference genome with HISAT2
echo "Indexing reference genome with HISAT2..."
hisat2-build -p 20 $HISAT2_INDEX_DIR/$REFERENCE_GENOME $HISAT2_INDEX_DIR/index

# HISAT2
echo "Running HISAT2..."
for base in $OUTPUT_DIR/trimmed/*_R1_paired.fastq; do
    base=$(basename $base "_R1_paired.fastq")
    hisat2 -p 15 -x $HISAT2_INDEX_DIR/index \
    -1 $OUTPUT_DIR/trimmed/${base}_R1_paired.fastq \
    -2 $OUTPUT_DIR/trimmed/${base}_R2_paired.fastq \
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
output_txt="$OUTPUT_DIR/prepDE/file_list.txt"  # Name of the file list that will be generated
> $output_txt  # Create or clear the text file

for bamfile in $OUTPUT_DIR/samtools/*_sorted.bam; do
    base=$(basename $bamfile "_sorted.bam")
    /path/to/stringtie/stringtie -e -B -p 16 -G $HISAT2_INDEX_DIR/$GFF3_FILE -o $OUTPUT_DIR/prepDE/${base}.gtf $bamfile
    echo "${base}.gtf $OUTPUT_DIR/prepDE/${base}.gtf" >> $output_txt  # Add the line to the text file
done

# prepDE.py
echo "Running prepDE.py..."
python $SCRIPTS_DIR/prepDE.py3 -i $output_txt -g $OUTPUT_DIR/prepDE/gene_count_matrix.csv -t $OUTPUT_DIR/prepDE/transcript_count_matrix.csv

echo "Pipeline finished!"
