#!/bin/bash

#---------------------------------------------#
#--------- Install SRA toolkits --------------#
#---------------------------------------------#

SRATOOLKIT_DIR=$(find . -maxdepth 1 -type d -name "sratoolkit.*" | head -n 1)

if [ -d "$SRATOOLKIT_DIR" ]; then
    echo ""
else
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xvzf sratoolkit.current-ubuntu64.tar.gz
    rm sratoolkit.current-ubuntu64.tar.gz
fi

#-----------------------------------------#
#--------- run SRA toolkits --------------#
#-----------------------------------------#

# Get the name of the extracted directory (should start with "sratoolkit.")
SRATOOLKIT_DIR=$(find . -maxdepth 1 -type d -name "sratoolkit.*" | head -n 1)

# Add toolkit binaries to PATH
export PATH="$PWD/$SRATOOLKIT_DIR/bin:$PATH"

# Check if fasterq-dump is now found
which fasterq-dump
#fasterq-dump --help

#------------------------------------------#
#--------- Install demo data --------------#
#------------------------------------------#

fasterq-dump SRR12345678 --split-files --outdir ./01_raw --threads 8

#Split the .sra files into:
#   SRR12345678_1.fastq <=> read 1
#   SRR12345678_2.fastq <=> read 2
#=> paired-end read

# Compress _1 files into .fq.gz to save space
parallel -j 4 --bar 'gzip {}; mv {}.gz ./01_raw/{/.}.fq.gz' ::: ./01_raw/*_1.fastq

# Compress _1 files into .fq.gz to save space
parallel -j 4 --bar 'gzip {}; mv {}.gz ./01_raw/{/.}.fq.gz' ::: ./01_raw/*_2.fastq



