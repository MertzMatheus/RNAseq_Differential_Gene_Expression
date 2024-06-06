# RNA-Seq Pipeline

This repository contains a bash script for processing RNA-Seq data using various bioinformatics tools. The pipeline includes quality control with FastQC, trimming with Trimmomatic, alignment with HISAT2, BAM file manipulation with Samtools, and expression quantification with StringTie. The script is designed to handle different sequencing technologies, including Illumina and Nanopore.

## Directory Structure

- **input_dir**: Directory containing the input FASTQ files.
- **output_dir**: Directory where the output files will be stored.
- **trimmomatic_dir**: Directory of Trimmomatic.
- **hisat2_index_dir**: Directory containing the HISAT2 indexed reference genome.
- **reference_genome**: Name of the reference genome file.
- **gff3_file**: Name of the GFF3 annotation file.
- **scripts_dir**: Directory containing the `prepDE.py` script.

## Requirements

Ensure the following tools are installed and available in your PATH:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [HISAT2](http://daehwankimlab.github.io/hisat2/)
- [Samtools](http://www.htslib.org/)
- [StringTie](https://ccb.jhu.edu/software/stringtie/)
- [Python3](https://www.python.org/)

### Installation

#### FastQC
Download and install FastQC from the [FastQC project page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

#### Trimmomatic
Download Trimmomatic from the [Trimmomatic page](http://www.usadellab.org/cms/?page=trimmomatic).

#### HISAT2
Download and install HISAT2 from the [HISAT2 page](http://daehwankimlab.github.io/hisat2/).

#### Samtools
Download and install Samtools from the [Samtools page](http://www.htslib.org/).

#### StringTie
Download and install StringTie from the [StringTie page](https://ccb.jhu.edu/software/stringtie/).

#### Python3
Ensure Python3 is installed on your system. Download it from the [Python official site](https://www.python.org/).

## How to Run

1. Clone this repository:
   ```bash
   git clone https://github.com/MertzMatheus/RNAseq_Pipeline.git
   cd RNAseq_Pipeline
Claro! Aqui está um exemplo detalhado de `README.md` para o seu diretório no GitHub, incluindo manuais dos pacotes utilizados e detalhamento de cada processo, além de mencionar diferentes opções de sequenciamento, como Illumina e Nanopore.

2. Modify the paths in the `RNAseq_Pipeline.sh` script as necessary for your setup.

3. Run the script:
   ```bash
   bash RNAseq_Pipeline.sh
   ```

The pipeline will process your RNA-Seq data and generate output files in the specified directories.

## Detailed Workflow

### 1. Quality Control with FastQC

FastQC performs quality control checks on raw sequence data. It provides a modular set of analyses that you can use to identify potential problems.

**Command:**
```
fastqc -t 15 -o $OUTPUT_DIR/fastqc $INPUT_DIR/*.fastq
```

### 2. Trimming with Trimmomatic

Trimmomatic trims adapter sequences and low-quality bases from the reads.

**Command:**
```
java -jar /path/to/Trimmomatic-0.39/trimmomatic.jar PE -threads 15 \
$sample $sample_R2 \
$OUTPUT_DIR/trimmed/${base}_R1_paired.fastq $OUTPUT_DIR/trimmed/${base}_R1_unpaired.fastq \
$OUTPUT_DIR/trimmed/${base}_R2_paired.fastq $OUTPUT_DIR/trimmed/${base}_R2_unpaired.fastq \
ILLUMINACLIP:/path/to/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:20 TOPHRED33
```

### 3. Indexing the Reference Genome with HISAT2

HISAT2 indexes the reference genome to prepare it for alignment.

**Command:**
```
hisat2-build -p 20 $HISAT2_INDEX_DIR/$REFERENCE_GENOME $HISAT2_INDEX_DIR/index
```

### 4. Alignment with HISAT2

HISAT2 aligns the RNA-Seq reads to the reference genome.

**Command:**
```
hisat2 -p 15 -x $HISAT2_INDEX_DIR/index \
-1 $OUTPUT_DIR/trimmed/${base}_R1_paired.fastq \
-2 $OUTPUT_DIR/trimmed/${base}_R2_paired.fastq \
-S $OUTPUT_DIR/hisat2/${base}.sam
```

### 5. BAM File Manipulation with Samtools

Samtools converts SAM files to BAM files, sorts them, and indexes them.

**Command:**
```
samtools view -bS $samfile > $OUTPUT_DIR/samtools/${base}.bam
samtools sort $OUTPUT_DIR/samtools/${base}.bam -o $OUTPUT_DIR/samtools/${base}_sorted.bam
samtools index $OUTPUT_DIR/samtools/${base}_sorted.bam
```

### 6. Quantification with StringTie

StringTie assembles the transcripts and estimates their abundance.

**Command:**
```
/path/to/stringtie/stringtie -e -B -p 16 -G $HISAT2_INDEX_DIR/$GFF3_FILE -o $OUTPUT_DIR/prepDE/${base}.gtf $bamfile
```

### 7. Generating Count Matrices with prepDE.py

The `prepDE.py` script generates gene and transcript count matrices.

**Command:**
```
python $SCRIPTS_DIR/prepDE.py3 -i $output_txt -g $OUTPUT_DIR/prepDE/gene_count_matrix.csv -t $OUTPUT_DIR/prepDE/transcript_count_matrix.csv
```

## Sequencing Technologies

### Illumina Sequencing

Illumina sequencing produces high-throughput short reads. This pipeline is optimized for processing paired-end Illumina data.

### Nanopore Sequencing

Nanopore sequencing produces long reads, which may require different preprocessing steps. Adjustments to the pipeline may be necessary to handle Nanopore data, such as using different trimming tools optimized for long reads.

## License

This project is licensed under the terms of the [MIT License](LICENSE).
```

Adapte os caminhos e outras configurações de acordo com suas necessidades específicas. Esse `README.md` fornece uma descrição detalhada dos processos e ferramentas utilizados no pipeline, além de orientações sobre como configurar e executar o script, abrangendo diferentes tecnologias de sequenciamento.
