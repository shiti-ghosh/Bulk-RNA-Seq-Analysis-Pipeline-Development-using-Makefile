# RNA-Seq Pipeline Directory Structure

    rna_seq_pipeline/
    │── README.md                  # Documentation for the pipeline
    │── dir_struct.pdf             # PDF representation of the directory structure
    │── env.yaml                   # Conda environment configuration file
    │── makefile                   # Makefile for automated pipeline execution
    │
    ├── data/                      # Data directory
    │   ├── input/                 # Raw input data
    │   │   ├── pair_end/          # Paired-end sequencing data
    │   │   │   ├── .gitkeep       # Placeholder to keep directory in Git
    │   │   ├── ref/               # Reference genome and annotation files
    │   │   │   ├── annotation/    # Annotation files
    │   │   │   ├── genome/grch38/ # Human genome reference (GRCh38)
    │   │   │   ├── transcriptome/ # Transcriptome reference
    │   ├── output/                # Processed output data
    │   │   ├── pair_end/          # Processed paired-end data
    │   │   │   ├── test/          # Test directory for output results
    │   │   │   │   ├── aligned/   # Aligned reads in BAM format
    │   │   │   │   ├── counts/    # Gene expression count data
    │   │   │   │   ├── fastqc/    # FastQC reports
    │   │   │   │   ├── multiqc/   # MultiQC reports
    │   │   │   │   ├── trimmed/   # Trimmed read outputs
    │   │   │   │   ├── .gitkeep   # Placeholder for version control
    │   ├── metadata.tsv           # Metadata file containing sample details
    │
    ├── logs/                      # Log files for debugging and tracking pipeline execution
    │   ├── .gitkeep               # Placeholder file to keep directory in Git
    │
    ├── script/                    # Directory for analysis scripts
    │   ├── .gitkeep               # Placeholder file
    │   ├── step1_fastqc.sh        # Step 1: Quality control with FastQC
    │   ├── step2_trimming.sh      # Step 2: Read trimming
    │   ├── step3_alignment.sh     # Step 3: Alignment using HISAT2
    │   ├── step4_sam_to_bam.sh    # Step 4: Convert SAM to BAM
    │   ├── step5_featurecounts.sh # Step 5: Gene expression quantification
    │   ├── step6_deseq2.R         # Step 6: Differential expression analysis
    │   ├── step7_visualization.R  # Step 7: Visualization (volcano plot, heatmap)
