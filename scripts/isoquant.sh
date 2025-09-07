conda activate isoquant
isoquant.py \
  --output ./results \
  --complete_genedb \
  --genedb ./reference/gencode.vM25.annotation.gtf \
  --reference ./reference/GRCm38.p6.genome.fa \
  --bam ./bam/OVA_TG_hifi.bam ./bam/PBS_TG_hifi.bam \
  --labels OVA_TG PBS_TG \
  --prefix TG.hifi \
  --data_type pacbio \
  --matching_strategy precise \
  --splice_correction_strategy default_pacbio \
  --model_construction_strategy default_pacbio \
  --transcript_quantification all \
  --gene_quantification all \
  --report_novel_unspliced true \
  --threads 200 \
  --check_canonical \
  --sqanti_output \
  --count_exons

