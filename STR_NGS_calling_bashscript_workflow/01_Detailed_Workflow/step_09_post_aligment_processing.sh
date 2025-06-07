#!/bin/bash

ref_dir="./00_ref"
align_dir="./03_align"
stats_dir="./04_bam_stats_flagstat"

export ref_dir align_dir stats_dir

files=( $(find "$align_dir" -type f -name "*.bam") )

post_align_process() {
    
    file="$1"
    base=$(basename "$file" .bam)

    # 1) Sort the BAM file
    #gatk SortSam \
    #    -I "$align_dir/${base}.bam" \
    #    -O "$align_dir/${base}_sorted.bam" \
    #    --SORT_ORDER coordinate
    samtools sort "$align_dir/${base}.bam" --threads 10 -o "$align_dir/${base}_sorted.bam"

    # 3) Mark Duplicates
    #gatk MarkDuplicates \
    #   --INPUT "$align_dir/${base}_sorted.bam" \
    #   --OUTPUT "$align_dir/${base}_dedup.bam" \
    #   --METRICS_FILE "$align_dir/${base}_MarkDup.metrics"

    # 4) Adding or Replacing Read Groups
    #gatk AddOrReplaceReadGroups \
    #    -I "$align_dir/${base}_dedup.bam" \
    #    -O "$align_dir/${base}_dedup_rg.bam" \
    #    --RGID "${base//_trimmed_align/}" \
    #    --RGLB lib1 \
    #    --RGPL ILLUMINA \
    #    --RGPU unit1 \
    #    --RGSM "${base//_trimmed_align/}"

    # 5) Index BAM file
    #gatk BuildBamIndex \
    #    -I "$align_dir/${base}_sorted.bam" \
    #    -O "$align_dir/${base}_sorted.bam.bai"
    #   -I "$align_dir/${base}_dedup_rg.bam" \
    #   -O "$align_dir/${base}_dedup_rg.bam.bai"
    samtools index -@ 10 "$align_dir/${base}_sorted.bam"

    # (optional) 6.1 Run BaseRecalibrator (creates a recalibration table)
    # gatk BaseRecalibrator \
    #    -I "$align_dir/${base}_output_with_rg.bam" \
    #    -R reference.fasta \
    #    --known-sites known_sites.vcf \
    #    -O recal_data.table

    # (optional) 6.2 Apply the recalibration to the BAM
    # gatk ApplyBQSR \
    #    -R reference.fasta \
    #    -I "$align_dir/${base}_output_with_rg.bam \
    #    --bqsr-recal-file recal_data.table \
    #    -O recal.bam
    
    # (optional) 7) samtools stats and flagstat for BAM file
    samtools stats --threads 10 \
        "$align_dir/${base}_sorted.bam" \
        > "$stats_dir/${base}_sorted_stats.txt"
    #    "$align_dir/${base}_dedup_rg.bam" \
    #    > "$stats_dir/${base}_dedup_rg_stats.txt"

    samtools flagstat --threads 10 \
        "$align_dir/${base}.bam" \
        > "$stats_dir/${base}_sorted_flagstat.txt"
    #     "$align_dir/${base}_dedup_rg.bam" \
    #     > "$stats_dir/${base}_dedup_rg_flagstat.txt"

    ## --threads 10 means using 10 threads for this task
}

export -f post_align_process

parallel -j 10 --bar post_align_process ::: "${files[@]}"