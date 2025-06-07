#!/bin/bash

# Define list of required directories
folders=(
    "./00_ref"
    "./01_raw"
    "./01_raw/qc_checked"
    "./02_trim"
    "./02_trim/qc_checked"
    "./03_align"
    "./04_bam_stats_flagstat"
    "./05_calling_annotation"
    "./05_calling_annotation/calling_dataframe"
)

# Loop through each folder and check/create
for folder in "${folders[@]}"; do
    if [ -d "$folder" ]; then
        echo "Directory already exists: $folder"
    else
        echo "Creating directory: $folder"
        mkdir -p "$folder"
    fi
done
