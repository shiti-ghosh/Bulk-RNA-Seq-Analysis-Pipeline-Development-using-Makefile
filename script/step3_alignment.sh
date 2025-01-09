#!/bin/bash
# Step 3: Alignment with HISAT2

READS_DIR="data/output/pair_end/test/trimmed"
ALIGN_DIR="data/output/pair_end/test/aligned"
INDEX_DIR="data/input/ref/genome/grch38"
THREADS=8

mkdir -p "$ALIGN_DIR"

echo "Aligning reads with HISAT2..."
for R1 in ${READS_DIR}/*_R1_val_1.fq.gz; do
    SAMPLE=$(basename $R1 _R1_val_1.fq.gz)
    R2=${READS_DIR}/${SAMPLE}_R2_val_2.fq.gz

    if [[ -f "$R2" ]]; then
        hisat2 -x ${INDEX_DIR}/genome \
            -1 $R1 -2 $R2 \
            -S ${ALIGN_DIR}/${SAMPLE}_aligned.sam
    fi
done

echo "Alignment completed!"

