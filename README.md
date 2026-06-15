# C. difficile Bulk RNA-seq Analysis Pipeline

This repository contains a Python-based pipeline for the quantification and downstream analysis of bulk RNA-seq data from *Clostridioides difficile* T6 strain.

---

## 1. Overview

(coming soon)

---

## 2. From FASTQ to BAM and basic read statistics

## 2.1 Build reference genome index

Before mapping RNA-seq reads, an index of the reference genome must be created. This allows the aligner to efficiently map sequencing reads to the genome.

We use the FASTA file (`.fna`) from the genome annotation set as the reference.

```bash
# Move into genome directory (optional)
cd Genomes

# Build bowtie2 index
bowtie2-build DMG2301622_Cdiff.fna Cdiff_index
```
## 2.2 Align reads to the reference genome (FASTQ → BAM)

After building the reference genome index, sequencing reads are aligned to the genome using `bowtie2`. The input consists of paired-end FASTQ files, and the output is a SAM file that is then converted into BAM format for efficient processing.

```bash
# Perform alignment using bowtie2
bowtie2 -x Genomes/Cdiff_index \
        -1 fastq/A10_METRO_E1_ZKRN260013229-1A_23H7NJLT3_L8_1.fq.gz \
        -2 fastq/A10_METRO_E1_ZKRN260013229-1A_23H7NJLT3_L8_2.fq.gz \
        | samtools view -bS - > mapped_reads.bam

# To get the total number of maps reads for normalization, run
samtools view -c -F 4 mapped_reads.bam

# Keep note of this number for later in the pipeline.

## 2.3 Sort and index aligned reads (BAM processing)

After alignment, the resulting BAM file is not yet ordered in genomic coordinate space. Sorting is required so that reads are arranged by their mapping position along the reference genome. This is necessary for efficient querying and visualization.

Once sorted, the BAM file is indexed to enable fast access to specific genomic regions.

```bash
# Sort the BAM file by genomic coordinates
samtools sort -o mapped_reads_sorted.bam mapped_reads.bam

# Index the sorted BAM file
samtools index mapped_reads_sorted.bam

# Once this is done, you can extract reads which map to specific loci (the example region corresponds to an rRNA operon)
samtools view -h -o mapped_reads_sorted_rRNA_operon.sam \
    mapped_reads_sorted.bam \
    "gnl|Prokka|CdiffT6_1:2971000-2977000"


