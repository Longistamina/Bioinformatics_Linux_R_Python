from Bio import SeqIO
import pandas as pd
import json
import re

dict_raw = json.load(open("ref_dir/raw_haplogroup_motif.json", "r"))
# {
#     "mt-MRCA": [
#         "16129A",
#         "16182M",
#         "16183M",
#         "16187T",
#         "16189C",
#         "16223T",
#         "16230G",
#          ......]

for record in SeqIO.parse("ref_dir/rcrs.fasta", "fasta"):
    rCRS = list(record.seq)

dict_standardized = dict()

'''This script only uppercases all motif and keeps the IUPAC nucleotide codes unchanged'''
'''It also converts all "del" into "gap" to facilitate following processing steps'''

for haplogroup, motif_list in dict_raw.items():
    dict_standardized[haplogroup] = dict()

    for motif in motif_list:
        pos, variant = re.search(r"(\d+\.?\d*)([a-zA-Z]+)", motif).groups()      
        dict_standardized[haplogroup][pos] = dict()

        if re.search(r"\d+\.{1}\d+", pos) is not None:
            ref = "gap"
        else:
            ref_index = int(pos) - 1
            ref = rCRS[ref_index]

        if variant == "del":
            variant = "gap" 
        else:
            variant = variant.upper()
        
        dict_standardized[haplogroup][pos].update({"ref": ref, "alt": variant})


with open("ref_dir/standardized_haplogroup_motif.json", "w") as file:
      json.dump(dict_standardized, file, indent=4)