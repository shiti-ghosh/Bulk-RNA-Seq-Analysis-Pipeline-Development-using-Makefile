# Variables
SCRIPTS_DIR = script
FASTQC_SCRIPT = $(SCRIPTS_DIR)/step1_fastqc.sh
TRIM_SCRIPT = $(SCRIPTS_DIR)/step2_trimming.sh
ALIGN_SCRIPT = $(SCRIPTS_DIR)/step3_alignment.sh
SAM_TO_BAM_SCRIPT = $(SCRIPTS_DIR)/step4_sam_to_bam.sh
FEATURECOUNTS_SCRIPT = $(SCRIPTS_DIR)/step5_featurecounts.sh
DESEQ2_SCRIPT = $(SCRIPTS_DIR)/step6_deseq2.R
VISUALIZATION_SCRIPT = $(SCRIPTS_DIR)/step7_visualization.R


# Targets
all: fastqc trimming alignment sam_to_bam featurecounts deseq2 visualization

fastqc:
	bash $(FASTQC_SCRIPT)

trimming: fastqc
	bash $(TRIM_SCRIPT)

alignment: trimming
	bash $(ALIGN_SCRIPT)

sam_to_bam: alignment
	bash $(SAM_TO_BAM_SCRIPT)

featurecounts: sam_to_bam
	bash $(FEATURECOUNTS_SCRIPT)
deseq2: featurecounts
	Rscript $(DESEQ2_SCRIPT)

visualization: deseq2
	Rscript $(VISUALIZATION_SCRIPT)
