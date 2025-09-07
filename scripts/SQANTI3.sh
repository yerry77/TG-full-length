#!/bin/bash
# SQANTI3 analysis pipeline for TG full-length transcriptome
# Note: Paths are placeholders. Replace them with your actual local paths.

# Activate conda environment
conda activate SQANTI3.env

# Set PYTHONPATH for cDNA_Cupcake utilities
export PYTHONPATH=$PYTHONPATH:/path/to/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/path/to/cDNA_Cupcake/

# Create output directories
mkdir -p ./SQANTI3_output

# Step 1: Run SQANTI3 QC
python /path/to/SQANTI3/sqanti3_qc.py \
    ./input/TG.transcript_models.gtf \   # corrected transcript GTF from IsoQuant
    ./reference/gencode.vM25.annotation.gtf \   # reference GTF
    ./reference/GRCm38.p6.genome.fa \           # reference genome FASTA
    --CAGE_peak ./reference/mouse.refTSS_v3.1.mm10.bed \
    --polyA_motif_list ./reference/mouse_and_human.polyA_motif.txt \
    -o TG \
    -d ./SQANTI3_output/QC_result/ \
    --cpus 16 \
    --report both \
    --short_reads ./input/short_reads.fofn

# Step 2: Run IsoAnnotLite
python /path/to/SQANTI3/utilities/IsoAnnotLite_SQ3.py \
    ./SQANTI3_output/QC_result/TG_corrected.gtf \
    ./SQANTI3_output/QC_result/TG_classification.txt \
    ./SQANTI3_output/QC_result/TG_junctions.txt \
    -o ./SQANTI3_output/IsoAnnotLite.gff3 \
    -novel -stdout ./SQANTI3_output

# Step 3: Filter isoforms with SQANTI3 rules
python /path/to/SQANTI3/sqanti3_filter.py rules \
    --isoAnnotGFF3 ./SQANTI3_output/IsoAnnotLite.gff3 \
    --isoforms ./SQANTI3_output/QC_result/TG_corrected.fasta \
    --gtf ./SQANTI3_output/QC_result/TG_corrected.gtf \
    --faa ./SQANTI3_output/QC_result/TG_corrected.faa \
    -o TG \
    -d ./SQANTI3_output/filtered_result \
    -j ./path/to/filter/filter_default.json \
    ./SQANTI3_output/QC_result/TG_classification.txt

# Step 4: Rescue filtered isoforms
python /path/to/SQANTI3/sqanti3_rescue.py rules \
    --isoforms ./SQANTI3_output/QC_result/TG_corrected.fasta \
    --gtf ./SQANTI3_output/filtered_result/TG.filtered.gtf \
    -g ./reference/gencode.vM25.annotation.gtf \
    -f ./reference/GRCm38.p6.genome.fa \
    -k ./SQANTI3_output/QC_result/TG_classification.txt \
    -e all \
    --mode full \
    -o TG \
    -d ./SQANTI3_output/rescue_result \
    -j ./path/to/filter/filter_default.json \
    ./SQANTI3_output/filtered_result/TG_RulesFilter_result_classification.txt

