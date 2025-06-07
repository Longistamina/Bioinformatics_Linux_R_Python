#!/bin/bash

#Download SRR32354557 sequence
wget https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR32354557

#Add gzip suffix (.gz) to use gzip -d to unzip
mv 'fastq?acc=SRR32354557' 'fastq?acc=SRR32354557.gz'
gzip -d 'fastq?acc=SRR32354557.gz'

#Add .fastq suffix
mv 'fastq?acc=SRR32354557' 'SRR32354557.fastq'
