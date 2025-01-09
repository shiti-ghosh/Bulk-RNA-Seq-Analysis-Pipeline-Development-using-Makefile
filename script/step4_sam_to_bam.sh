#!/bin/bash
# Step 4: Convert SAM to BAM and sort

ALIGN_DIR="data/output/pair_end/test/aligned"

echo "Converting SAM to BAM and sorting..."
for SAM_FILE in ${ALIGN_DIR}/*.sam; do
    SAMPLE=$(basename "$SAM_FILE" "_aligned.sam")
    
    # Generate sorted BAM filename
    BAM_FILE=${ALIGN_DIR}/${SAMPLE}_aligned.sorted.bam
    
    # Convert SAM to BAM and sort
    samtools view -Sb "$SAM_FILE" | samtools sort -o "$BAM_FILE"
    
    # Optional: Remove the original SAM file
    rm "$SAM_FILE"
    
    echo "Processed $SAM_FILE -> $BAM_FILE"
done

echo "SAM to BAM conversion and sorting completed!"

