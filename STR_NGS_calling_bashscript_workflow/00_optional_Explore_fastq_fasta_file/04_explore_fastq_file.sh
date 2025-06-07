#!/bin/bash

#----------------------------------------------------------#
#---------------- for .fq.gz compressed file --------------#
#----------------------------------------------------------#

# zcat SRR32354557.fq.gz | awk 'END{print NR, NF}'

echo " "

# This command prints the number of fields (NF) and the number of lines (NR) in the FASTQ file
awk 'END{print NR, NF}' SRR32354557.fastq
#=> Output: 2701256 1

#--------------------------------------------#
#--------------------------------------------#
echo " "

# This command counts the number of sequences (lines starting with '@') in the FASTQ file
echo " "
awk '/^@/ {count++} END {print count}' SRR32354557.fastq
#The /..../ is the regular expression
#So the /^@/ is used to match the line starting with "@"
#=> Output: 679184

#--------------------------------------------#
#--------------------------------------------#
echo " "

# Print sequence ID and length
awk 'NR % 4 == 1 {print $1, $3}' SRR32354557.fastq | head -n 20

#As explained in file 03_explain_fastq_fastq_file.sh, a fastq sequence consists of 4 lines
#Use NR % 4 == 1 to get each sequence
#    NR is number of record, so if NR = 1, meaning the ID line of the first sequence, then if NR % 4 = 1 % 4 = 1
#    get this line to access $1 and $3
#    $1 is field 1 is the ID, 
#    $3 is field 3 is the length

#If NR = 5, then NR % 4 = 5 % 4 = 1
#   then get this NR line as the ID line of sequence 2

#Do the same for other equence
echo " "

#Or can use /^@/: !awk '/^@/ {print $1, $3}' SRR32354557.fastq | head -n 20
awk '/^@/ {print $1, $3}' SRR32354557.fastq | head -n 20

#=> Output: @SRR32354557.1 length=151
#           @SRR32354557.1 length=150
#           @SRR32354557.2 length=151
#           @SRR32354557.2 length=151

#--------------------------------------------#
#--------------------------------------------#
echo " "

# Print the ID and length of sequences that are shorter than 40
awk 'NR % 4 == 1 {id = $1} NR % 4 == 2 {if (length($0) < 40) print id, "length="length($0)}' SRR32354557.fastq | head

# NR % 4 == 2 is the 2nd line in each sequenche block, it is the sequence line ATGCCT....
# $0 means everything, meaning the whole sequence, so length($0) to get the sequence length

#=> Output: @SRR32354557.854 length=36
#	    @SRR32354557.854 length=36
#	    @SRR32354557.2650 length=35
#	    @SRR32354557.2650 length=35

#--------------------------------------------#
#--------------------------------------------#
echo " "

# Find the longest sequence in FASTQ
awk 'BEGIN {len_max = 0} NR % 4 == 1 {id = $1} NR % 4 == 2 {if (len_max < length($0)) {len_max=length($0); len_max_id=id; seq=$0}} END {print "Longest_sequence: "len_max_id, "\n"seq, "\nlength="len_max}' SRR32354557.fastq
#=> Output:
# Longest_sequence: @SRR32354557.1 
# GGAGTTAAAATTGCACTTGCATGCAAAATCGCTTCTGTATGTGTATCATCAAAATCAAGTGCATGAGAAGCAATGCCATCTAGCATGGCTGCGTATATAGGATCTAATTTTTGATCATTTATCCAAATAAGAGTATTTTTTGTAAGGGGTA 
# length=151

echo " "

# Find the shortest sequence in FASTQ
awk 'NR == 2 {len_min = length($0)} NR % 4 == 1 {id = $1} NR % 4 == 2 {if (len_min > length($0)) {len_min=length($0); len_min_id=id; seq=$0}} END {print "Shortest_sequence: "len_min_id, "\n"seq, "\nlength="len_min}' SRR32354557.fastq
#=> Output:
# Shortest_sequence: @SRR32354557.2650 
# GGGTGGAGTCATTGAAATTTCGCATTTAAGACAGT 
# length=35

#--------------------------------------------#
#--------------------------------------------#
echo " "

# Count A nucleotides in each sequence of SRR32354557.fastq
awk 'NR % 4 == 1 {id = $1} NR % 4 == 2 {for (i=1; i<=length($0); i++) {if (substr($0, i, 1)=="A") A++}; print id, "| Total A = "A}' SRR32354557.fastq | head

# substr($0, i, 1) means extract a substring start with i-th character with length = 1,  from the entire $0 input line
# For example, if $0 = ATTTCGAA, and i=2
# Then substr($0, i, 1) = T

#=> Output: @SRR32354557.1 | Total A = 48
#           @SRR32354557.1 | Total A = 106

#--------------------------------------------#
#--------------------------------------------#
echo " "

# Count A, T, G, C nucleotides in FASTQ
awk 'NR % 4 == 1 {id = $1} NR % 4 == 2 {for (i=1; i<=length($0); i++) {if (substr($0, i, 1)=="A") A++;  else if (substr($0, i, 1)=="T") T++; else if (substr($0, i, 1)=="G") G++; else C++}; print id, "| Total A = "A, "| Total T = "T, "| Total G = "G, "| Total C = "C}' SRR32354557.fastq | head

#=> Output: @SRR32354557.1 | Total A = 48 | Total T = 50 | Total G = 31 | Total C = 22
#           @SRR32354557.1 | Total A = 106 | Total T = 109 | Total G = 43 | Total C = 43

#-------------------------------------------#
#-------------------------------------------#
echo " "

# Remove '@' from sequence headers
awk '/^@/ {gsub("@", ""); print $0}' SRR32354557.fastq | head

#/^@/ to find the records starting with "@"
#Then {gsub("@", ""); print $0} to replace "@" with "", and print out all content of that record
# The ; is delimiter between commands inside the {}
#head to print out only first 5 records

#=> Output: SRR32354557.1 1 length=151
#	    SRR32354557.1 1 length=150
#	    SRR32354557.2 2 length=151
#           SRR32354557.2 2 length=151

#--------------------------------------------#
#--------------------------------------------#
echo " "

# Find sequences containing ATTTAAT
awk 'NR % 4 == 1 {id = $1} NR % 4 == 2 && /ATTTAAT/ {print id, "seq:"$0}' SRR32354557.fastq | head

#The /..../ is the regular expression
#So, /ATTAAT/ is used to match the sequence containing the pattern "ATTTAAT"

#=> Output: @SRR32354557.4 seq:GGTAAAAAACACTTTTGGGATATAAAATTTAAAATCAAAGCTCGTGTTTATCCGTCTTTTAAACAGCATTTAATCACTTTTTTTGGGACAAAACGGATATTTTGGGATATTTGGGATTTTTTACACTTTAAATTAGG
#	    @SRR32354557.5 seq:TAGTTACATCAATAAGAAAACATAAGAATACTTTAAGCAAAAATGTAACTATCCTAGAGGGTAAAAAACACTTTTGGGATATAAAATTTAAAATCAAAGCTCGTGTTTATCCGTCTTTTAAACAGCATTTAATCACT

#--------------------------------------------#
#--------------------------------------------#
echo " "

# Convert FASTQ to FASTA
awk 'NR % 4 == 1 {print ">"substr($1,2)} NR % 4 == 2 {print $0}' SRR32354557.fastq | head

#substr($1, 2): The substr() function returns a substring of the input string. In this case:
#    + The first argument is $1, which is the entire first field of the line (usually the header, starting with a @ in FASTQ files).
#    + The second argument is 2, which means "start from the second character."
#=> This effectively removes the leading @ character from the sequence header. So, if a line starts like @ABC123, substr($1, 2) would return ABC123.

#print ">"substr(...) to add the ">" at the begining of header line, representing the FASTA file

#=> Output: 
# >SRR32354557.1
# GGAGTTAAAATTGCACTTGCATGCAAAATCGCTTCTGTATGTGTATCATCAAAATCAAGTGCATGAGAAGCAATGCCATCTAGCATGGCTGCGTATATAGGATCTAA
# >SRR32354557.1
# TCCTTTAACTTTTCAATCACTTCTTGCATTTGCATTTTAAATCCTTTGATTTTAACTTCAAAGAATTTAAACATAAAATTTCAAAAAGCAATTAAAAACTTTTGATT
