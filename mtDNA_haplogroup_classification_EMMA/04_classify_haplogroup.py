import math
import numpy as np
import json
import os
import logging
from argparse import ArgumentParser
import shutil
from multiprocessing import Pool

logger = logging.getLogger(__name__)

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
    

    def load_samples_profile (self, json_input_file):
        '''
        Extract each sample from the .json file.
        Then, convert all its variant profile into SAM format.
        And then cast them into a dataframe.
        '''
        sample_profile_path = os.path.join(self.input_dir, json_input_file)
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

        
    def classify_haplogroup(self, sample):
        haplogroup_result = dict()
        haplogroup_result[sample] = dict()

        for haplogroup in self.haplogroup_motif.keys():
            cost_sum = 0

            haplogroup_pos_set = set(self.haplogroup_motif[haplogroup].keys())
            sample_pos_set = set(self.dict_samples_profile[sample].keys())

            intersect_pos = haplogroup_pos_set & sample_pos_set # mustation that exists in both sample and haplogroup
            #diff_sample_pos = sample_pos_set - haplogroup_pos_set # mutation that sample has but haplogroup does not

            for pos in self.dict_samples_profile[sample].keys(): # mustation that exists in both sample and haplogroup

                if (pos == "ranges") or (str(math.floor(float(pos))) in Disregard_InDels):
                    continue

                sample_base = self.dict_samples_profile[sample][pos].get('alt')

                if pos in intersect_pos:
                    motif_base = self.haplogroup_motif[haplogroup][pos].get('alt')
                else:
                    motif_base = self.dict_samples_profile[sample][pos].get('ref')
                    
                if (set(IUPAC_CODE.get(motif_base)) >= set(IUPAC_CODE.get(sample_base))): # check for heteroplasmy in motif
                    continue                      
                elif (self.cost_value.get(pos) is not None):
                    cost_sum += self.assign_existed_cost(pos, sample_base, motif_base)
                else:
                    cost_sum += HaplogroupClassifier.assign_new_cost(sample_base, motif_base)

            haplogroup_result[sample][haplogroup] = cost_sum
        
        haplogroup_result[sample] = dict(sorted(haplogroup_result[sample].items(), key = lambda x: x[1])) # sort the dict ascending

        return haplogroup_result
        

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
    args = parse_args()
    classifier = HaplogroupClassifier(
        input_dir = args.input, 
        output_dir = args.output, 
        ref_dir = args.reference, 
        cores = args.cores
    )
    
    try:
        classifier.load_haplogroup_motif()
        classifier.load_cost_value()

        for json_file in os.listdir(classifier.input_dir):
            sample_name = json_file.replace(".json", "")
            classifier.dict_samples_profile[sample_name] = dict()
            classifier.dict_samples_profile[sample_name].update(classifier.load_samples_profile(json_file))

        #classifier.classify_haplogroup("242749")
        ###### classify ######
        arg_list = [(sample,) for sample in classifier.dict_samples_profile.keys()]
        
        with Pool(processes = classifier.cores) as pool:
            result_list = pool.starmap(func = classifier.classify_haplogroup, iterable = arg_list)
        
        haplogroup_result_dict = dict()
        for sample_haplogroup in result_list:
            haplogroup_result_dict.update(sample_haplogroup)
        
        for sample, haplgroup_data in haplogroup_result_dict.items():
            json_path = os.path.join(classifier.output_dir, f"{sample}.json")
            with open(json_path, "w") as json_file:
                json.dump(haplgroup_data, json_file, indent = 4)
        ###### classify ######
        
        logger.info("All samples processed successfully") 
    
    except Exception as e:
        logger.error(f"Error in main {e}")

if __name__ == '__main__':
    main()
