Long-read RNA-seq maps the isoform landscape of the trigeminal ganglion in allergic rhinitis

## Abstract
This repository contains scripts, processed data, and analysis pipelines for investigating novel isoforms in the trigeminal ganglion using full-length transcriptome sequencing.

This study identified novel isoforms, quantified transcript expression levels, and provided analysis pipelines and plotting scripts.

## Introduction
The trigeminal ganglion (TG) plays a crucial role in allergic rhinitis (AR) by innervating the nasal mucosa through the branched projection of its sensory neurons.  It mediates AR-associated symptoms such as nasal itching and sneezing, likely through TRP channel activation and altered electrophysiological properties. However, the underlying sensory transcriptomes remain incompletely understood. While conventional short-read RNA-sequencing (RNA-seq) has been widely used for transcriptome analysis, its short read lengths limit accurate characterization of complex transcripts, including novel isoforms and long non-coding RNAs. In contrast, long-read full-length transcriptome sequencing provides complete transcript structures, enabling the discovery of novel genes and splice variants, but suffers from lower quantitative accuracy. To overcome these limitations, we employed an integrative approach, combing short-read and long-read sequencing technologies to achieve a comprehensive and precise analysis of the TG transcriptome in the context of AR-like mouse model.

## Data Description
Data Source

Organism: Mus musculus

Tissue: Trigeminal ganglion (TG)

Experimental model: Allergic rhinitis (AR) mouse model and control group

Sequencing technology: PacBio HiFi (Iso-Seq, full-length transcriptome sequencing)

## Workflow
For a detailed view of the analysis workflow, please refer to the following PDF document:

[Download Analysis Workflow PDF](https://github.com/yerry77/TG-full-length/blob/main/analysis%20workflow.pdf)
## Results 
novel isoforms gtf file.[novel isoforms.gtf](https://github.com/yerry77/TG-full-length/blob/main/novel%20isoforms.gtf)  

TPM values for novel isoforms calculated using RSEM. [novel isoform TPM use RSEM.tsv](https://github.com/yerry77/TG-full-length/blob/main/novel%20isoform%20TPM%20use%20RSEM.tsv)

Figure1 data [Figure1 data.xlsx](https://github.com/yerry77/TG-full-length/blob/main/Figure1%20data.xlsx)

