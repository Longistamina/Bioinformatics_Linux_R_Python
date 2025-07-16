import math
import numpy as np
import pandas as pd
import json
import os
from loguru import logger
from argparse import ArgumentParser
import shutil
import re
from multiprocessing import Pool

IUPAC_CODE = {
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "T", "G", "C"),
    "gap": ("gap",)
}

VARIANT_SET = {"A", "C", "G", "T", "gap"}

Disregard_InDels = ("16193", "309", "455", "463", "573", "960", "5899", "8276", "8285")

class HaplogroupClassifier:

    def __init__ (self, input_dir, output_dir, ref_dir, cores = 2):
        '''
        input_dir: path to input directory (should content input in .json file)
        output_dir: path to output directory storing haplogroup classification results
        ref_path: path to reference .json file storing information of haplogroups in SAM format
                                                                                     (SAM = String-based Alignment for MtDNA)
        cores: number of CPU cores to use for parallel processing
        '''
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.ref_dir = ref_dir
        self.cores = cores
        self.dict_samples_profile = dict()

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)


    def load_haplogroup_motif (self):
        haplogroup_motif_path = os.path.join(self.ref_dir, "standardized_haplogroup_motif.json")
        with open(haplogroup_motif_path, "r") as json_file:
            self.haplogroup_motif = json.load(json_file)
        
        return self.haplogroup_motif
    

    def load_cost_value (self):
        cost_value_path = os.path.join(self.ref_dir, "cost_value.json")
        with open(cost_value_path, "r") as json_file:
            self.cost_value = json.load(json_file)
        
        return self.cost_value
    

    def _load_sample_profile (self, json_sample_file):
        '''
        Extract each sample from the .json file.
        Then, convert all its variant profile into SAM format.
        And then cast them into a dataframe.
        '''
        sample_profile_path = os.path.join(self.input_dir, json_sample_file)
        with open(sample_profile_path, "r") as json_file:
            dict_sample_profile = json.load(json_file)
        
        return dict_sample_profile
    

    @staticmethod
    def assign_new_cost(sample_base, motif_base):
        cost_value = 0

        if set(IUPAC_CODE.get(motif_base)) >= set(IUPAC_CODE.get(sample_base)):
            return cost_value

        else:
            sample_base_tuple = IUPAC_CODE.get(sample_base)
            for motif_base in IUPAC_CODE.get(motif_base):
                fr_dict = dict()

                for fluc_base in list(VARIANT_SET - {motif_base}):
                    transition_check = ({motif_base, fluc_base} == {"A", "G"}) or ({motif_base, fluc_base} == {"C", "T"})
                    
                    match transition_check:
                        case True:
                            fluctuation_rate = 10**(-6)
                        case _:
                            fluctuation_rate = 10**(-9)
                    
                    fr_dict[fluc_base] = fluctuation_rate

                stable = 1 - sum(fr_dict.values())
                
                for sample_base in sample_base_tuple:
                    if sample_base == motif_base:
                        continue
                    cost_value += np.log10(stable / fr_dict[sample_base])/3

        return cost_value


    def assign_existed_cost(self, pos, sample_base, motif_base):
        sample_base_tuple = IUPAC_CODE.get(sample_base) # handle heteroplasmy of sample base
        motif_base_tuple = IUPAC_CODE.get(motif_base) # handle heteroplasmy of motif base

        cost_value = 0
        for sample_base in sample_base_tuple:
            for motif_base in motif_base_tuple:
                if sample_base == motif_base:
                    continue
                else:
                    cost_value += self.cost_value[pos][motif_base].get(sample_base)
        return cost_value

        
    def classify_haplogroup(self, json_sample_file):
        haplogroup_result = dict()
        dict_sample_profile = self._load_sample_profile(json_sample_file)

        for haplogroup in self.haplogroup_motif.keys():
            cost_sum = 0

            haplogroup_pos_set = set(self.haplogroup_motif[haplogroup].keys())
            sample_pos_set = set(dict_sample_profile.keys())

            intersect_pos = haplogroup_pos_set & sample_pos_set # mustation that exists in both sample and haplogroup
            #diff_sample_pos = sample_pos_set - haplogroup_pos_set # mutation that sample has but haplogroup does not

            for pos in dict_sample_profile.keys(): # mustation that exists in both sample and haplogroup

                if (pos == "ranges") or (str(math.floor(float(pos))) in Disregard_InDels):
                    continue

                sample_base = dict_sample_profile[pos].get('alt')

                if pos in intersect_pos:
                    motif_base = self.haplogroup_motif[haplogroup][pos].get('alt')
                else:
                    motif_base = dict_sample_profile[pos].get('ref')
                    
                if (set(IUPAC_CODE.get(motif_base)) >= set(IUPAC_CODE.get(sample_base))): # check for heteroplasmy in motif
                    continue                      
                elif (self.cost_value.get(pos) is not None):
                    cost_sum += self.assign_existed_cost(pos, sample_base, motif_base)
                else:
                    cost_sum += HaplogroupClassifier.assign_new_cost(sample_base, motif_base)

            haplogroup_result[haplogroup] = cost_sum
        
        haplogroup_result = dict(sorted(haplogroup_result.items(), key = lambda x: x[1])) # sort the dict ascending

        classify_path = os.path.join(self.output_dir, json_sample_file)
        with open(classify_path, "w") as json_file:
            json.dump(haplogroup_result, json_file, indent = 4)
        

#-------------------------------------------------------#                
#----------------- main() function ---------------------#
#-------------------------------------------------------#

def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(description="Generate JSON files from variant data")
    parser.add_argument("-i", "--input", required=True, help="Path to input directory containing sample files (input_dir/standardized)")
    parser.add_argument("-o", "--output", required=True, help="Path to output directory (classify_dir/)")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference directory containing cost_value.json file (ref_dir/)")
    parser.add_argument("-c", "--cores", type=int, default=2, required=False, help="The number of CPU cores to be used for parallel computing, default is 2. (Example: set -c 8 means using 8 cores to process 8 different samples at a time)")
    return parser.parse_args()

def main():
    # classifier = HaplogroupClassifier(
    #     input_dir = "input_dir/standardized", 
    #     output_dir = "classify_dir", 
    #     ref_dir = "ref_dir/", 
    #     cores = 14
    # )
    try:
        args = parse_args()
        classifier = HaplogroupClassifier(
            input_dir = args.input, 
            output_dir = args.output, 
            ref_dir = args.reference, 
            cores = args.cores
        )
        classifier.load_haplogroup_motif()
        classifier.load_cost_value()

        #classifier.classify_haplogroup("242749")
        
        ###### classify ######
        arg_list = [(json_sample_file,) for json_sample_file in os.listdir(classifier.input_dir)]
        
        with Pool(processes = classifier.cores) as pool:
            pool.starmap(func = classifier.classify_haplogroup, iterable = arg_list)
        ###### classify ######
        
        logger.info("All samples processed successfully") 
    
    except Exception as e:
        logger.error(e)

if __name__ == '__main__':
    main()


# python3 python3 04_classify_haplogroup.py -i input_dir/standardized -o classify_dir/ -r ref_dir/ -c 10