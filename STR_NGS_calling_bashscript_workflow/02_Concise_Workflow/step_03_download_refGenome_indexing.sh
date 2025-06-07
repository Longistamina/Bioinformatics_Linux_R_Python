#!/bin/bash

# If already have human reference genome, skip this file and steps

#---------------------------------------------------#
#------- download reference genome -----------------#
#---------------------------------------------------#

# Download hg38
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -O ./00_ref/hg38.fa.gz
gunzip -q ./00_ref/hg38.fa.gz


#---------------------------------------------------#
#------- indexing reference genome -----------------#
#---------------------------------------------------#

# Use bwa-mem2
bwa-mem2 index ./00_ref/hg38.fa 
# After running this command, BWA will generate index files (.amb, .ann, .bwt, .pac, .sa) in the same folder as the reference FASTA.

# Use samtools
samtools faidx ./00_ref/hg38.fa #resulted in hg38.fa.fai


