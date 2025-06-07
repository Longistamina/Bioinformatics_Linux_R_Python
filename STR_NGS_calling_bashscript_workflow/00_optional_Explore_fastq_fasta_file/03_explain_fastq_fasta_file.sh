#---------------------------------------------#
#------------- fastq file --------------------#
#---------------------------------------------#

FASTQ File

Purpose: Stores raw sequencing reads from platforms like Illumina.
Contents: Contains both the sequence and quality scores for each base.

Format: 4 lines per read
###
@sequence_ID
AGTCTAGCTAGCTAGCTAGC
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**
###

Line 1: @ + sequence identifier
Line 2: nucleotide sequence
Line 3: + (optionally followed by the same sequence ID)
Line 4: ASCII-encoded quality scores (same length as sequence)

Used for:
	Raw data from sequencing machines
	Input to tools like FastQC, Trimmomatic, BWA, STAR
	
Example:
###
@SRR123456.1 1 length=76
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAAGGGTGCCCGATAGGACCTCGGATGAAGTG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SRR123456.2 2 length=76
CGTAAGCTCTGAACTTCAGGGTCAAGGAAGGAACTGTGATGTGTGGCTGAATAGTGGTTCTATTTC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
...
###


#---------------------------------------------#
#------------- fasta file --------------------#
#---------------------------------------------#

FASTA File

Purpose: Stores nucleotide or amino acid sequences.
Contents: Just the sequence, optionally with a description.

Format:
###
>sequence_ID description (optional)
AGTCTAGCTAGCTAGCTAGC
AGCTAGCTAGCTA
###

Used for:
	Reference genomes
	Protein sequences
	Input for aligners (e.g., BWA, BLAST)
	
Example:
###
>NC_012920.1 Homo sapiens mitochondrion, complete genome
GATCACAGGTCTATCACCCTATTAAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCT
GGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGTCAAGGTGTAGCCCATGAGGCCAAATATCA
TTCTGAGGGGAGACCCGATCTATTTTACCACAGACAAACCTACGCCAAAATCCATTTCACTATCAG
...
###
