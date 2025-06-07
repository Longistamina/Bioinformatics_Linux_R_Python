#------------------------------------------------------#
#----------------- install STRsensor ------------------#
#------------------------------------------------------#

# cd wanted_directory

# git init
# git clone https://github.com/ChenHuaLab/STRsensor.git

# chmod +x ./STRsensor/STRsensor_1.2.2_x64-Linux

# mv ./STRsensor/STRsensor_1.2.2_x64-Linux ./STRsensor/str-sensor

# export PATH=$PATH:./STRsensor

# str-sensor #check if command is found

#-------------------------------------------------------#

ref_dir="./00_ref"
align_dir="./03_align"
calling_dir="./05_calling_annotation"

for file in "$align_dir"/*_aligned_sorted.bam; do
    echo "$file" >> "$calling_dir/list_bam_file.txt"
done

str-sensor -i "$calling_dir/list_bam_file.txt" \
          -f "$ref_dir/hg38.fa" \
          -r "$ref_dir/STR_hg38.txt" \
          -n 10 \
          -t 20 \
          -m 20 \
          -o ./05_calling_annotation

rm "$calling_dir/list_bam_file.txt"

python3 ./generate_report.py


# Usage: STRsensor -i <bamlist.txt> -r <region.txt> -f <genome.fa> -o <output_dir> [options]
# Options:
#        -h|--help             print help infomation

# [Required]
#        -i|--infile     FILE  bam file list, one sample per line [.txt]
#        -r|--region     FILE  region file of the STR locus [.txt]
#        -f|--fasta      FILE  reference genome (must be same as the mapping process used) [.fa]
#        -o|--outpath    PATH  the path for output file

# [Optional input parameters]
#        -s|--stutter    FILE  the stutter parameters file for each STR locus [.txt]
#        -q|--frequency  FILE  the frequency file for each STR locus on each potential allele [.txt]

# [Optional output parameters]
#        -p|--params_out       output parameters learned by given samples (including
#                              stutter and allele frequency)

# [Optional filter rules]
#        -c|--min_prob   FLOAT the minimal probability allowed to call an allele [0.00]
#        -n|--min_reads  INT   the minimum reads needed to genotype a STR locus for an individual [10]
#        -m|--mis_match  INT   the maximum mismatch(mismatch/insert/delete) bases allowed in both
#                              3' and 5' flanking sequence [2]

# [Optional]
#        -a|--allow_dup        duplication is allowed in target and amplicon sequencing, but not in WGS
#        -t|--threads    INT   the number of threads used [1]


# ************** NOTICE **************
# USE SHORT FLAGS LIKE -i -r f instead of long flags --infile --region ...
# If the number of reads of a seuqence is TOO LOW (< 50 000 reads) => can be CORE DUMPED error
