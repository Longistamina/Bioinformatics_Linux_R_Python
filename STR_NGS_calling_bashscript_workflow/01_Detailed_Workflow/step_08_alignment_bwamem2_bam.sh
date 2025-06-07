#!/bin/bash

ref="./00_ref/hg38.fa"
trim_dir="./02_trim"
align_dir="./03_align"

#Export align_dir (and other needed variables) so that GNU parallel can access them:
export ref trim_dir align_dir

files=( $(find "$trim_dir" -type f -name "*.fq.gz") )

# Each sample consumes 25 threads (-t 25)
# Each time processes 10 samples in parallel (-j 4)
# => consume 100 threads (100 cpu cores)

align_to_bam() {
    file="$1"

    # Skip if this is the second read of a paired-end file to avoid duplicate processing
    if [[ "$file" == *_2.fq.gz ]]; then
        echo "Skipping $file because it is the second read of a pair."
        return
    fi

    # For paired-end, check if the paired file exists
    if [[ "$file" == *_1.fq.gz ]] && [[ -f "${file/_1.fq.gz/_2.fq.gz}" ]]; then
        echo "Detected paired-end reads for $file. Running paired-end alignment..."
        r1="$file"
        r2="${r1/_1.fq.gz/_2.fq.gz}"
        base=$(basename "$r1" _1.fq.gz)
        bwa-mem2 mem -t 6 "$ref" "$r1" "$r2" | samtools view -Sb - > "$align_dir/${base}_aligned.bam"
    else
        echo "Detected single-end reads for $file. Running single-end alignment..."
        base=$(basename "$file" .fq.gz)
        bwa-mem2 mem -t 6 "$ref" "$file" | samtools view -Sb - > "$align_dir/${base}_aligned.bam"
    fi
}

export -f align_to_bam

parallel -j 4 --bar align_to_bam ::: "${files[@]}"
