#------------------------------------------------------------#
#----------------- Awk and Gawk introduction ----------------#
#------------------------------------------------------------#

# awk is a powerful tool in bash script for text processing task
awk -W version #check awk version, and check if awk is available

# gawk is a GoogleColab version of awk
sudo apt-get update
sudo apt-get install gawk #install gawk

gawk --verion

# awk and gawk syntax
awk 'criteria {action}' input_file.txt #input_file can be a .txt, .fastq or .fasta file
gawk 'criteria {action}' input_file.txt

# awk and gawk with pipeline
echo "My name is Long" | awk 'print $4' #output: Long
echo "My name is Long" | gawk 'print $4' #output: Long


#------------------------------------------------------#
#----------------- -F: field separator ----------------#
#------------------------------------------------------#

echo "My name is Long" | awk '{print $4}' 
# By default, awk understands " " (space character) as field separator
# $4 means field number 4th
#=> Output: Long

echo "My,name,is,Long" | awk '{print $4}' 
# Not thing is printed out since field separator is no longer " "
#=> Output: "" (nothing)

echo "My,name,is,Long" | awk -F "," '{print $4}' 
# -F "," to define "," as field separator
#=> Output: Long

echo "My,name,is,Long" | awk -F "," '{print $2, $4}'
# Print out only the $2 and the $4 field
#=> Output: name Long


#--------------------------------------------------------------------------#
#----------------- BEGIN action _ main action _ END action ----------------#
#--------------------------------------------------------------------------#

#'BEGIN{action_beg} {action_main} END{action_end}' input_file
# to execute action_beg, then action_main with input_file, then action_end

echo -e "Apple\nBanana\nCherry" | awk 'BEGIN{print "Fruits list:"} {print $1} END{print "End of list."}'

# echo -e => enable special Escape Sequences like \n = newline, \t = horizontal tab, \v = vertical tab, ...

# BEGIN{print "Fruits list:"} => begin action, print the string "Fruits list:" before the main action (don't need input)
# {print $1}                  => main action, take the input from echo command and print out the $1 field
# END{print "End of list."}   => end action, print the string "End of list." after the execution of main action

#=> Output: Fruits list:
#           Apple
#           Banana
#           Cherry
#           End of list.


#-------------------------------------------------------#
#----------------- .awk file and awk -f ----------------#
#-------------------------------------------------------#

# A .awk file contains awk commands that can be executed via awk -f command
# awk -f => to execute awk commands in a .awk file

echo -e "Hello Word\nAWK is great" > file.txt
# Generate text input in file.txt

echo -e 'BEGIN {print "Header"} {print $1, $2} END {print "Footer"}' > script.awk
# Generate awk commands into script.awk file

# cat file.txt
# cat script.awk

awk -f script.awk practice_data/file.txt
# Excute awk commands in script.awk, with input is file.txt
# awk -f Read the AWK program source from the file program-file,  instead of  from  the  first command line argument.

#=> Output: Header
#           Hello Word
#           AWK is
#           Footer


#---------------------------------------------------#
#----------------- NR, NF, $NR, $NF ----------------#
#---------------------------------------------------#

# NR is the number of records (rows)
# NF is the number of fields (columns)

# $NF to access the last field (column)

awk 'END{print NR, NF}' practice_data/awk_example_data.tsv
# Print out the number of records and fields of the input file
#=> Output: 38 9

awk '{print NR, NF}' practice_data/awk_example_data.tsv
#=> Output: 1 9
#           2 9
#           3 9
#           ...
#           38 9

awk '{print $NF}' practice_data/awk_example_data.tsv
# Print out the last field of the input file
#=> Output: Gene
#           WASH7P
#           TUBB8
#           WBP1LP10
#           ....
#           MIR1302-9HG

awk 'END{print}' practice_data/awk_example_data.tsv
awk 'END{print $0}' practice_data/awk_example_data.tsv
# Both print out the contents of the last record of the input file
# $0 means print out everything
#=> Output: chr9	27657	30891	25.0	23.4	0.004	3	0.001	MIR1302-9HG

awk 'END{print $NR, $NF}' practice_data/awk_example_data.tsv
# Print out the last record of the last field
#=> Output: MIR1302-9HG


#----------------------------------------------------------------------#
#----------------- Use criteria to query data with awk ----------------#
#----------------------------------------------------------------------#

awk 'criteria {action}' input_file.txt #input_file can be a .txt, .fastq or .fasta file

# awk comparison operators:
# Operator	Description	         	Example
# ==	        Equal to	         	$1 == "apple"
# !=	        Not equal to	         	$2 != 0
# <	        Less than		 	$3 < 100
# >	        Greater than		        $3 > 50
# <=	        Less than or equal to	        $1 <= 10
# >=		Greater than or equal to	$2 >= 5
# ~	        Matches regex pattern		$1 ~ /apple/
# !~	        Does not match regex pattern	$2 !~ /^b/


# awk logical operators:
# Operator	Description	         	Example
# &&	        AND (true if all true)         	awk '($1 > 10) && ($2 == "pass") { print $0 }' data.txt
# ||	        OR (false if all false)        	awk '($3 == "high") || ($4 < 50) { print $0 }' data.txt
# !	        Inverts the truth	 	awk '!($5 ~ /error/) { print $0 }' data.txt


awk '$1=="chr12" {print $0}' practice_data/awk_example_data.tsv
# Print out the rows/records of awk_example_data.tsv file where the field $1 equal to "chr12"
#=> Output: chr12	12310	13501	21.4	19.8	0.01	4	0.0032	DDX11L8
#           chr12	14522	32015	29.0	27.3	0.21	1	0.02	WASH8P
#           chr12	36602	38133	24.5	22.0	0.001	5	0.001	FAM138D


echo "This is the table with q-values <= 0.001" | gawk '{print $0}'
awk 'BEGIN{getline; print $0} ($8<=0.001){print $0} END{print "All files are processed!"}' practice_data/awk_example_data.tsv | column -t
# {getline; print $0} to get the 1st line of the input file (the column line), and print out all columns of this line
# condition is ($8 <= 0.001)
# column -t => to parse the output into column format with tab character (-t) as separator


#-----------------------------------------------------------#
#----------------- Awk execute many actions ----------------#
#-----------------------------------------------------------#

awk '{action 1; action 2; action 3; ....}' input_file #To execute multiple actions, each action in a separate line

awk -F "," '$1=="Yes" {yes++} $1=="No" {no++} END {print "People with HeartDisease: " yes; print "People without HeartDisease: " no}' practice_data/heart_2020_cleaned.csv

# if $1=="Yes" then execute action {yes++} or yes = yes + 1
# if $1=="No" then execute action {no++} or no = no + 1
# then {print "People with HeartDisease: " yes; print "People without HeartDisease: " no}
#=> Output:
#   People with HeartDisease: 27373
#   People without HeartDisease: 292422


#-----------------------------------------------------------#
#----------------- gsub for substituting text --------------#
#-----------------------------------------------------------#

gsub(REGEXP, REPLACEMENT, TARGET)
# If TARGET is omitted, gsub operates on $0 (the whole line).

echo "The cat chased another cat." | awk '{gsub(/cat/, "dog"); print $0}'
#=> Output: The dog chased another dog.

echo "ID1234Name5678" | awk '{gsub(/[0-9]/, ""); print $0}'
#=> Output: IDName

echo "hello_world_this_is_awk" | awk '{gsub(/_/, " "); print $0}'
#=> Output: hello world this is awk

echo "apple, banana, cherry" | awk -F ", " '{gsub(/a/, "A", $2); print $0}'
#=> Output: apple bAnAnA cherry


#---------------------------------------------------------------------#
#----------------- substr for getting substring of text --------------#
#---------------------------------------------------------------------#

substr(string, start, length)
# string: the input string.
# start: the position (1-based) where to start extracting.
# length: (optional) number of characters to extract. If omitted, extracts to the end of the string.

echo -e "apple\nbanana\ncherry" | awk '{ print substr($0, 1, 3) }'
#=> Output: app
#	    ban
#	    che

echo "gene1234" | awk '{ print substr($1, 2, 4) }'
#=> Output: ene1

echo "chromosome" | awk '{ print substr($1, 6) }'
#=> Output: some
