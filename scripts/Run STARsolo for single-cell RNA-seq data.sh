#!/bin/bash
# Script: Run STARsolo for single-cell RNA-seq data.sh
# Description: Run STARsolo for single-cell RNA-seq BAM file to obtaining splice junction information from single-cell data

# User-configurable variables
THREADS=80
GENOME_DIR="/path/to/STAR_index" #STAR index data of the reference genome
INPUT_BAM="/path/to/possorted_genome_bam.bam" # BAM file of the single-cell RNA-seq data
WHITELIST="/path/to/737K-august-2016.txt" # Barcode whitelist file
OUTPUT_DIR="./SRR11925964" # Output directory

# Create output directory
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}" || exit 1

# Run STARsolo
/data/p/bin/STAR \
  --runThreadN ${THREADS} \
  --genomeDir ${GENOME_DIR} \
  --soloType CB_UMI_Simple \
  --readFilesIn ${INPUT_BAM} \
  --readFilesCommand "samtools view -F 0x100" \
  --readFilesType SAM SE \
  --soloInputSAMattrBarcodeSeq CR UR \
  --soloInputSAMattrBarcodeQual CY UY \
  --soloCBwhitelist ${WHITELIST} \
  --soloFeatures Gene SJ \
  --soloOutFileNames SRR11925964

echo "STARsolo analysis completed. Output saved to: ${OUTPUT_DIR}"

# The above provides a sample running script. You can run other samples by changing the sample.
