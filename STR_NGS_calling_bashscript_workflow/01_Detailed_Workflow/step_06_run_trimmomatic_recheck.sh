#!/bin/bash

# find ~/ -name "NexteraPE-PE.fa"
#Copy the path: /home/longdp/anaconda3/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa

#------------------------------------------------------------------#
#------------------- Trimming with trimmomatic --------------------#
#------------------------------------------------------------------#

files_raw=( $(find ./01_raw -type f -name "*.fq.gz") )

trimming() {

  file="$1"
  base=$(basename "$file" .fq.gz)

  # Paired-end case
  if [[ "$file" == *_1.fq.gz ]]; then
    r1="$file"
    r2="${file/_1.fq.gz/_2.fq.gz}"
    pair_base=$(basename "$r1" _1.fq.gz)

    if [[ -f "$r2" ]]; then
      echo "=> Trimming paired-end: $pair_base"

      trimmomatic PE \
        -threads 25 \
        "$r1" "$r2" \
        ./02_trim/"$pair_base"_trimmed_1.fq.gz ./02_trim/"$pair_base"_unpaired_1.fq.gz \
        ./02_trim/"$pair_base"_trimmed_2.fq.gz ./02_trim/"$pair_base"_unpaired_2.fq.gz \
        SLIDINGWINDOW:4:10 \
        AVGQUAL:10 \
        MINLEN:20 \
        ILLUMINACLIP:/home/longdp/anaconda3/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10:8:3:true
        #The unpaired files store the reads that was discarded during trimming
        #=> this helps ensure all paired-end reads still have their mate
    else
      echo "!!! Missing pair for $r1, skipping."
    fi

  # Single-end case (and ignore already-processed _2.fastq)
  elif [[ "$file" != *_2.fq.gz ]]; then
    echo "=> Trimming single-end: $base"

    trimmomatic SE \
      -threads 25 \
      "$file" ./02_trim/"$base"_trimmed.fq.gz \
      SLIDINGWINDOW:4:10 \
      AVGQUAL:10 \
      MINLEN:20 \
      ILLUMINACLIP:/home/longdp/anaconda3/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10:8:3:true
  fi
}

export -f trimming

parallel -j 4 --bar trimming ::: "${files_raw[@]}"
  

#-------------------------------------------------------------#
#---------------- recheck QC after trimming ------------------#
#-------------------------------------------------------------#

files_trimmed=( $(find ./02_trim -type f -name "*.fq.gz") )

qc_after_trimming() {

    file="$1"

    # If it's paired-end (_1.fq.gz), handle both files
    if [[ "$file" == *_1.fq.gz ]]; then
        
        r1="$file"
        r2="${file/_1.fq.gz/_2.fq.gz}"
        pair_base=$(basename "$r1" _1.fq.gz)

        if [[ -f "$r2" ]]; then
            echo "=> Paired-end processing: $pair_base"
            fastqc "$r1" "$r2" -t 10 -o ./02_trim/qc_checked
        else
            echo "!!! Missing pair for $r1, skipping."
        fi

    elif [[ "$file" != *_2.fq.gz ]]; then
        
        # Skip R2s if they've already been handled
        base=$(basename "$file" .fq.gz)
        echo "=> Single-end processing: $base"
        fastqc "$file" -t 10 -o ./02_trim/qc_checked
    
    fi
}

export -f qc_after_trimming

parallel -j 10 --bar qc_after_trimming ::: "${files_trimmed[@]}"