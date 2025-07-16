import math
import json
import os
import re
from loguru import logger
from argparse import ArgumentParser
import shutil
from multiprocessing import Pool

Disregard_InDels = ("16193", "309", "455", "463", "573", "960", "5899", "8276", "8285")

class VariantHarmoniser ():
    
    def __init__(self, sample_dir, classify_dir, harmonise_dir, motif_path, cores = 2):
        self.sample_dir = sample_dir        # Input 1: variant profile of each sample
        self.classify_dir = classify_dir    # Input 2: haplogroup classification result of each sample
        self.motif_path = motif_path        # Input 3: haplogroup motif profile
        self.harmonise_dir = harmonise_dir  # Output: harmonised profile of each sample
        self.cores = cores # The number of CPU cores used for multiprocessing
        
        if os.path.exists(harmonise_dir):
            shutil.rmtree(harmonise_dir)
        os.makedirs(harmonise_dir)


    def load_haplogroup_motif (self):
        with open(self.motif_path, "r") as json_file:
            self.haplogroup_motif = json.load(json_file)
        
        return self.haplogroup_motif
    
    
    def _load_sample_profile (self, json_sample_file):
        sample_profile_path = os.path.join(self.sample_dir, json_sample_file)
        with open(sample_profile_path, "r") as json_file:
            dict_sample_profile = json.load(json_file)
            dict_sample_profile.pop('ranges')
        
        return dict_sample_profile
    
    
    def _extract_sample_haplogroup (self, json_sample_file):
        sample_haplogroup_path = os.path.join(self.classify_dir, json_sample_file)
        with open(sample_haplogroup_path, "r") as json_file:
            dict_sample_haplogroup = json.load(json_file)

        best_haplogroup = next(iter(dict_sample_haplogroup.keys())) # Get the first key of the dictionary, as well as the best haplogroup

        return best_haplogroup
    

    def harmonise_variant(self, json_sample_file):
        sample_profile = self._load_sample_profile(json_sample_file)
        best_haplogroup = self._extract_sample_haplogroup(json_sample_file)
        best_motif = self.haplogroup_motif[best_haplogroup]

        motif_pos_list = best_motif.keys()
        sample_pos_list = sample_profile.keys()

        intersect_pos = set(map(lambda pos: str(math.floor(float(pos))), motif_pos_list)) & set(map(lambda pos: str(math.floor(float(pos))), sample_pos_list))

        insert_pattern = r"\d+\.\d+"
        
        for motif_pos in motif_pos_list:
            motif_pos_floor = str(math.floor(float(motif_pos)))
 
            if motif_pos_floor in Disregard_InDels:
                continue
            
            elif (motif_pos_floor in intersect_pos) and (re.search(insert_pattern, motif_pos) is None):
                sample_profile[motif_pos]['alt'] = best_motif[motif_pos]['alt']

            elif re.search(insert_pattern, motif_pos):
                for sample_pos in sample_pos_list:

                    if sample_pos == motif_pos:
                        sample_profile[sample_pos]['alt'] = best_motif[motif_pos]['alt']
                    elif re.search(fr"{motif_pos_floor}\.\d+", sample_pos) and (float(sample_pos) > float(motif_pos)):
                        sample_profile.pop(sample_pos)
                    else:
                        continue
        
        harmonise_path = os.path.join(self.harmonise_dir, json_sample_file)
        with open(harmonise_path, "w") as json_file:
            json.dump(sample_profile, json_file, indent = 4)


#---------------------------------------------------#
#--------------- main() function -------------------#
#---------------------------------------------------#

def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(description="Generate JSON files from variant data")
    parser.add_argument("-s", "--sample", required=True, help="Path to sample directory containing sample files (input_dir/standardized)")
    parser.add_argument("-cls", "--classify", required=True, help="Path to classify directory (classify_dir/)")
    parser.add_argument("-m", "--motif", required=True, help="Path to standardized haplogroup motif file (ref_dir/standardized_haplogroup_motif.json)")
    parser.add_argument("-hmn", "--harmonise", required=True, help="Path to output harmonise directory (harmonise_dir/)")
    parser.add_argument("-c", "--cores", type=int, default=2, required=False, help="The number of CPU cores to be used for parallel computing, default is 2. (Example: set -c 8 means using 8 cores to process 8 different samples at a time)")
    return parser.parse_args()


def main():
    # harmoniser = VariantHarmoniser(
    #     sample_dir = "input_dir/standardized", 
    #     classify_dir = "classify_dir/", 
    #     motif_path = "ref_dir/standardized_haplogroup_motif.json", 
    #     harmonise_dir = "harmonise_dir/",
    #     cores = 2
    # )
    try:
        args = parse_args()
        harmoniser = VariantHarmoniser(
            sample_dir = args.sample, 
            classify_dir = args.classify, 
            motif_path = args.motif, 
            harmonise_dir = args.harmonise,
            cores = args.cores
        )
        harmoniser.load_haplogroup_motif()

        arg_list = [(json_sample_file,) for json_sample_file in os.listdir(harmoniser.sample_dir)]

        with Pool(processes = harmoniser.cores) as pool:
            pool.starmap(func = harmoniser.harmonise_variant, iterable = arg_list)
    
    except Exception as e:
        logger.error(e)


if __name__ == "__main__":
    main()


# python3 05_harmonise_variant.py -s input_dir/standardized -cls classify_dir/ -m ref_dir/standardized_haplogroup_motif.json -hmn harmonise_dir/ -c 10