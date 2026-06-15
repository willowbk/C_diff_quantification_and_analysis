# C. difficile Bulk RNA-seq Analysis Pipeline

This repository contains a Python-based pipeline for the quantification and downstream analysis of bulk RNA-seq data from *Clostridioides difficile*.

---

## 1. Overview

(coming soon)

---

## 2. From FASTQ to BAM and basic read statistics

### 2.1 Generate a BAM file from FASTQ

First, align reads from FASTQ files to a reference genome using a standard aligner such as `bowtie2` or `bwa`.

Example using `bowtie2`:

```bash
bowtie2 -x reference_genome_index \
        -1 sample_R1.fastq.gz \
        -2 sample_R2.fastq.gz \
        -S mapped_reads.sam




