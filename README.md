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

# To get the total number of mapped reads for normalization, run
samtools view -c -F 4 mapped_reads.bam

# Keep note of this number for later in the pipeline.
```
## 2.3 Sort and index aligned reads (BAM processing)

Sorting is required so that reads are arranged by their mapping position along the reference genome.

Once sorted, the BAM file is indexed to enable fast access to specific genomic regions.

```bash
# Sort the BAM file by genomic coordinates
samtools sort -o mapped_reads_sorted.bam mapped_reads.bam

# Index the sorted BAM file
samtools index mapped_reads_sorted.bam

# At this point, you can also generate a wiggle file in bigWig format to view in a genomebrowser,

bamCoverage \
    -b mapped_reads_sorted.bam \
    -o coverage.bw \
    --normalizeUsing CPM
    --binSize 1

# or for separate forward / reverse tracks:

bamCoverage \
    -b mapped_reads_sorted.bam \
    -o coverage_forward.bw \
    --normalizeUsing CPM \
    --filterRNAstrand forward \
    --binSize 1

bamCoverage \
    -b mapped_reads_sorted.bam \
    -o coverage_reverse.bw \
    --normalizeUsing CPM \
    --filterRNAstrand reverse \
    --binSize 1

# For WIG instead (typically larger), simply use

bamCoverage \
    -b mapped_reads_sorted.bam \
    -o coverage_forward.wig \
    --outFileFormat wig \
    --normalizeUsing CPM \
    --filterRNAstrand forward \
    --binSize 1

bamCoverage \
    -b mapped_reads_sorted.bam \
    -o coverage_reverse.wig \
    --outFileFormat wig \
    --normalizeUsing CPM \
    --filterRNAstrand reverse \
    --binSize 1

# Once this is done, you can extract reads which map to specific loci (the example region corresponds to an rRNA operon)
samtools view -h -o mapped_reads_sorted_rRNA_operon.sam \
    mapped_reads_sorted.bam \
    "gnl|Prokka|CdiffT6_1:2971000-2977000"

# Then move this final SAM file to the SAM_files folder
mv mapped_reads_sorted_rRNA_operon.sam SAM_files/

# It is also recommended to delete the previously created BAM files, as these can take up a lot of space.

```

## 2.4 Normalizing and plotting the mapped reads

Finally, the mapped reads can be visualized with the following Python script.

The parameters for which region to visualize are hard-coded at the start. 
These parameters are currently set for the example rRNA operon region.

```bash
python Plot_read_mappings.py

