#!/bin/bash
# Step 5: Generate counts using FeatureCounts

ALIGN_DIR="data/output/pair_end/test/aligned"
OUTPUT_COUNTS_DIR="data/output/pair_end/test/counts"
REF_ANNOTATION="data/input/ref/annotation/Homo_sapiens.GRCh38.113.gtf"

mkdir -p "$OUTPUT_COUNTS_DIR"

echo "Running FeatureCounts..."
featureCounts -p -a "$REF_ANNOTATION" -o "$OUTPUT_COUNTS_DIR/counts_summary.txt" ${ALIGN_DIR}/*_sorted.bam

echo "FeatureCounts completed!"

