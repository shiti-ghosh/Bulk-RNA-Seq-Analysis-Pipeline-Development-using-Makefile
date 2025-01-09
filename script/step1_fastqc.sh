#!/bin/bash
# Step 1: Run FastQC on raw data

INPUT_DIR="data/input/pair_end"
FASTQC_DIR="data/output/pair_end/test/fastqc/raw"
MULTIQC_DIR="data/output/pair_end/test/multiqc/raw"
THREADS=8

mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR"

echo "Running FastQC on raw data..."
while IFS=$'\t' read -r sample group subject; do
    fastq_file_r1="${INPUT_DIR}/${sample}_R1.fastq.gz"
    fastq_file_r2="${INPUT_DIR}/${sample}_R2.fastq.gz"

    if [[ -f "$fastq_file_r1" && -f "$fastq_file_r2" ]]; then
        fastqc -t "$THREADS" "$fastq_file_r1" "$fastq_file_r2" -o "$FASTQC_DIR"
    fi
done < "data/metadata.tsv"

echo "FastQC completed!"

echo "Running MultiQC on raw data..."
multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR"
