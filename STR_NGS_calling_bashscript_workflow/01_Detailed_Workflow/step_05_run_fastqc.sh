#!/bin/bash

#for file in `ls ./01_raw/*.fastq`; do fastqc $file -o ./01_raw/qc_checked/; done

files=( $(find ./01_raw -type f -name "*.fq.gz") )

qc_before_trimming() {

    file="$1"

    # If it's paired-end (_1.fq.gz), handle both files
    if [[ "$file" == *_1.fq.gz ]]; then
        
        r1="$file"
        r2="${file/_1.fq.gz/_2.fq.gz}"
        pair_base=$(basename "$r1" _1.fq.gz)

        if [[ -f "$r2" ]]; then
            echo "=> Paired-end processing: $pair_base"
            fastqc "$r1" "$r2" -t 10 -o ./01_raw/qc_checked
        else
            echo "!!! Missing pair for $r1, skipping."
        fi

    elif [[ "$file" != *_2.fq.gz ]]; then
        
        # Skip R2s if they've already been handled
        base=$(basename "$file" .fq.gz)
        echo "=> Single-end processing: $base"
        fastqc "$file" -t 10 -o ./01_raw/qc_checked
    
    fi
}

export -f qc_before_trimming

parallel -j 10 --bar qc_before_trimming ::: "${files[@]}"