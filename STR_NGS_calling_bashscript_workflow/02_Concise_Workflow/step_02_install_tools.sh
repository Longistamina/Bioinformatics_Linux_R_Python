#!/bin/bash

#Install conda based on this link: https://www.youtube.com/watch?v=7-naqq9fvZE&t=328s

#Activate conda virtual invironment



## fastqc is a tool for quality control of fastq file
## trimmomatic is a tool for trimming low-quality adapters or short sequences, help enhance reads quality
## fastp is a modern tool for QC and trimming, it combines all features of fastqc and trimmomatic
## bwa-mem2 is an alignment and index tool for SHORT READ, an imporved version of bwa
## minimap2 is like bwamem2 but for LONG READ
## samtools help convert SAM to BAM, sorting, indexing, do stats or flagstat
## gatk4 helps in post-alignment processing like mark duplicates, add readgroups, etc...

#-----------------------------------------------#
#---------- install all in one -----------------#
#-----------------------------------------------#

conda install -c  bioconda fastqc trimmomatic fastp minimap2 bwa-mem2 samtools gatk4

