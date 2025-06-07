#!/bin/bash

#Install conda based on this link: https://www.youtube.com/watch?v=7-naqq9fvZE&t=328s

#Activate conda virtual invironment

#-------------------------------------------#
#---------- install fastqc -----------------#
#-------------------------------------------#

# fastqc is a tool for quality control of fastq file

conda install -c bioconda fastqc

fastqc #check if installed successfully

#------------------------------------------------#
#---------- install trimmomatic -----------------#
#------------------------------------------------#

# trimmomatic is a tool for trimming low-quality adapters or short sequences
# trimmomatic helps enhance the quality of the sequence data before further analysis

conda install -c bioconda trimmomatic

trimmomatic #check if installed successfully

#------------------------------------------#
#---------- install fastp -----------------#
#------------------------------------------#

# fastp is a modern tool for QC and trimming, it combines all features of fastqc and trimmomatic

conda install -c bioconda fastp

fastp #check if installed successfully

#---------------------------------------------#
#---------- install bwa-mem2 -----------------#
#---------------------------------------------#

# bwa-mem2 is for SHORT READ
# an imporved version of bwa

conda install -c bioconda bwa-mem2
bwa-mem2

#---------------------------------------------#
#---------- install minimap2 -----------------#
#---------------------------------------------#

# minimap2 is for LONG READ
# minimap2 is a tool for sequences mapping and alignments like bwa or bwa-mem

# Note: minimap2 has replaced BWA-MEM for PacBio and Nanopore read alignment. 
# It retains all major BWA-MEM features, but is ~50 times as fast, more versatile, more accurate and produces better base-level alignment. 
# BWA-MEM2 is 50-100% faster than BWA-MEM and outputs identical alignments.

# Minimap2 is optimized for x86-64 CPUs

conda install -c bioconda minimap2

minimap2 #check if installed successfully


#---------------------------------------------#
#---------- install samtools -----------------#
#---------------------------------------------#

# The aligment step results in a .sam file

# samtools help convert SAM to BAM
# BAM is a binary compressed of SAM file, which is less heavier than SAM.
# Also helps in mapped reads processing

conda install -c bioconda samtools

samtools #check if installed successfully

#-----------------------------------------#
#---------- install gatk -----------------#
#-----------------------------------------#

# gatk helps in mark read duplications

conda install -c bioconda gatk4

gatk #check if installed successfully

#-----------------------------------------------#
#---------- install all in one -----------------#
#-----------------------------------------------#

# conda install -c  bioconda fastqc trimmomatic fastp minimap2 bwa-mem2 samtools gatk4

