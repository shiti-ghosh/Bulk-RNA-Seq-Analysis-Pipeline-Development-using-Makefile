# RNA-Seq Analysis Pipeline using Makefile

## Overview
This repository provides an RNA-Seq analysis pipeline for paired-end sequencing data of the human genome. 
The pipeline consists of multiple steps, from quality control to differential gene expression analysis, and is organized 
using a `Makefile` for efficient execution. Users can run all steps sequentially using `make` or execute individual steps as needed.
Begin by cloning the repository, if needed.

You may run the pipeline end to end or may choose steps and run them individually to have a better understanding of your outputs from each step.
For example, you may want to look at the QC reports from the first two steps to understand if "trimming" of sequences needed or not.

This pipeline can currently only take in paired end reads from Homo sapiens. \

Remember to edit the `metadata.tsv`  based on experiment conditions.

## Steps in the Pipeline

### Step 1: Create Virtual Environment and Install Dependencies
#### Significance: Ensures a controlled environment with all required dependencies.
```bash
conda env create -f env.yaml
conda activate rna_seq_pipeline
```

### Step 2: Run FastQC for Quality Control
#### Significance: Checks the quality of raw sequencing reads before processing.
```bash
make fastqc
```

### Step 3: Run MultiQC for Aggregated Quality Control Reports
#### Significance: Aggregates multiple FastQC reports into a single summary.
```bash
make multiqc
```

### Step 4: Trim Reads (Optional but Recommended)
#### Significance: Removes adapter sequences and low-quality bases to improve alignment accuracy.
```bash
make trimming
```

### Step 5: Align Reads with HISAT2
#### Significance: Maps trimmed reads to the human reference genome.
```bash
make alignment
```

### Step 6: Convert SAM to BAM and Sort
#### Significance: Converts alignment files to a more efficient BAM format and sorts them for downstream analysis.
```bash
make sam_to_bam
```

### Step 7: Generate Read Counts with FeatureCounts
#### Significance: Quantifies gene expression by counting mapped reads per gene.
```bash
make featurecounts
```

### Step 8: Perform Differential Gene Expression Analysis with DESeq2
#### Significance: Identifies significantly differentially expressed genes between conditions.
```bash
make deseq2
```

### Step 9: Visualization - Volcano Plot and Heatmap
#### Significance: Helps interpret results by highlighting differentially expressed genes.
```bash
make visualization
```

## Running the Entire Pipeline
To run all steps together using the Makefile:

```bash
make all
```

## Makefile Structure
Makefile program contains the pipeline steps in order. May be executed as a whole if you want to run the program without looking at the ouputs from the
intermediate steps.

## Running with Docker/Container
You can also containerize the pipeline using Docker:
```bash
docker build -t rnaseq_pipeline .
docker run -v $(pwd):/data rnaseq_pipeline
```

## Notes
- Ensure paths in the scripts match your directory structure.
- Modify `Makefile` to customize the workflow based on organism-specific requirements.
- If you need to change the reference genome, modify files in the `data/input/ref/` directory.
- The Conda environment is specified in `env.yaml` and includes all dependencies.

