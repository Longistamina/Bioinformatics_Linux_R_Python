#!/bin/bash

# fastp is a more modern tool
# fastp combines both fastqc and trimmomatic functions

#------------------------------------#
#---------- Install fastp -----------#
#------------------------------------#

#Config channel
## conda config --add channels defaults
## conda config --add channels bioconda
## conda config --add channels conda-forge

#Install fastp with bioconda
## conda install -c bioconda fastp

#Example for paired-end reads:
## fastp -i input_R1.fastq -I input_R2.fastq -o output_R1.fastq -O output_R2.fastq

#Example for single-end reads:
## fastp -i input.fastq -o output.fastq


#------------------------------------------------------#
#---------- quality control BEFORE trimming -----------#
#------------------------------------------------------#

if [ -d "/tmp/fastp_qc_dummy" ]; then
  echo "Directory /tmp/fastp_qc_dummy already exists, skipping creation."
else
  mkdir -p /tmp/fastp_qc_dummy
fi

# Each sample consumes 10 threads (--thread 10)
# Each time processes 10 samples in parallel (-j 10)
# => consume 100 threads (50 cpu cores)

files=( $(find ./01_raw -type f -name "*.fq.gz") ) #find all .fq.gz file in ./01_raw and store the names into list files


#Create qc_trimming() function
qc_before_trimming() {
  
  file="$1"

  # If it's paired-end (_1.fq.gz), handle both files
  if [[ "$file" == *_1.fq.gz ]]; then
    r1="$file"
    r2="${file/_1.fq.gz/_2.fq.gz}"
    pair_base=$(basename "$r1" _1.fq.gz)

    if [[ -f "$r2" ]]; then
      echo "=> Paired-end processing: $pair_base"
      fastp \
        -i "$r1" -I "$r2" \
        -o /tmp/fastp_qc_dummy/"$pair_base"_R1_dummy.fq.gz \
        -O /tmp/fastp_qc_dummy/"$pair_base"_R2_dummy.fq.gz \
        --disable_adapter_trimming \
        --thread 5 \
        --html ./01_raw/qc_checked/"$pair_base"_qc.html \
        --json ./01_raw/qc_checked/"$pair_base"_qc.json
    else
      echo "!!! Missing pair for $r1, skipping."
    fi

  elif [[ "$file" != *_2.fq.gz ]]; then
    # Skip R2s if they've already been handled
    base=$(basename "$file" .fq.gz)
    echo "=> Single-end processing: $base"
    fastp \
      -i "$file" \
      -o /tmp/fastp_qc_dummy/"$base"_dummy.fq.gz \
      --disable_adapter_trimming \
      --thread 5 \
      --html ./01_raw/qc_checked/"$base"_qc.html \
      --json ./01_raw/qc_checked/"$base"_qc.json
  fi
}

#Export to subshell
export -f qc_before_trimming

# Then run in parallel
parallel -j 4 --bar qc_before_trimming ::: "${files[@]}"


#-----------------------------------------------------#
#---------- quality control AFTER trimming -----------#
#-----------------------------------------------------#

# Each sample consumes 10 threads (--thread 10)
# Each time processes 10 samples in parallel (-j 10)
# => consume 100 threads (50 cpu cores)

#Create trimming_recheck() function
trimming_recheck() {
  file="$1"

  # If it's paired-end (_1.fq.gz), handle both files
  if [[ "$file" == *_1.fq.gz ]]; then
    r1="$file"
    r2="${file/_1.fq.gz/_2.fq.gz}"
    pair_base=$(basename "$r1" _1.fq.gz)

    if [[ -f "$r2" ]]; then
      echo "=> Paired-end processing: $pair_base"
      fastp \
        -i "$r1" -I "$r2" \
        -o ./02_trim/"$pair_base"_trimmed_1.fq.gz \
        -O ./02_trim/"$pair_base"_trimmed_2.fq.gz \
        --qualified_quality_phred 10 \
        --unqualified_percent_limit 50 \
        --length_required 20 \
        --n_base_limit 4 \
        --thread 5 \
        --html ./02_trim/qc_checked/"$pair_base"_trim_qc.html \
        --json ./02_trim/qc_checked/"$pair_base"_trim_qc.json
    else
      echo "!!! Missing pair for $r1, skipping."
    fi

  elif [[ "$file" != *_2.fq.gz ]]; then
    # Skip R2s if they've already been handled
    base=$(basename "$file" .fq.gz)
    echo "=> Single-end processing: $base"
    fastp \
      -i "$file" \
      -o ./02_trim/"$base"_trimmed.fq.gz \
      --qualified_quality_phred 10 \
      --unqualified_percent_limit 50 \
      --length_required 20 \
      --n_base_limit 4 \
      --thread 5 \
      --html ./02_trim/qc_checked/"$base"_trim_qc.html \
      --json ./02_trim/qc_checked/"$base"_trim_qc.json
  fi
}

#Export to subshell
export -f trimming_recheck

# Then run in parallel
parallel -j 4 --bar trimming_recheck ::: "${files[@]}"
