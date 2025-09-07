# TG Full-length Transcriptome Project
This project is used to analyze trigeminal ganglion (TG) full-length transcriptome sequencing data, including data processing scripts, result visualization and analysis workflow.
## use isoquant to generate long-read model
Use ultra to align the raw file of long-read sequencing to the reference genome to obtain a bam file, which serves as the input file for isoquant.
## use SQANTI3 to identify isoform and quality control
The output file of isoquant is used as the input file of SQANTI3 to perform quality control on the transcripts of the long-read model and identify novel isoforms.
