#!/bin/bash
# Step 2: Trimming with Trim Galore

INPUT_DIR="data/input/pair_end"
TRIM_DIR="data/output/pair_end/test/trimmed"
FASTQC_TRIM="data/output/pair_end/test/fastqc/trimmed"
MULTIQC_TRIM="data/output/pair_end/test/multiqc/trimmed"
THREADS=8

mkdir -p "$TRIM_DIR" "$FASTQC_TRIM" "$MULTIQC_TRIM"

echo "Trimming reads with Trim Galore..."
while IFS=$'\t' read -r sample group subject; do
    r1_in="${INPUT_DIR}/${sample}_R1.fastq.gz"
    r2_in="${INPUT_DIR}/${sample}_R2.fastq.gz"

    if [[ -f "$r1_in" && -f "$r2_in" ]]; then
        trim_galore --paired --cores "$THREADS" \
            --output_dir "$TRIM_DIR" \
            "$r1_in" "$r2_in"
    fi
done < "data/metadata.tsv"

echo "Trimming completed!"
echo "Running FASTQC on trimmed data..."
fastqc -t "$THREADS" "$TRIM_DIR"/*.fq.gz -o "$FASTQC_TRIM"

echo "Running MultiQC on trimmed data..."
multiqc "$FASTQC_TRIM" -o "$MULTIQC_TRIM"
