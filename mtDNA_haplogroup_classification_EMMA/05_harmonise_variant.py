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

        motif_pos = best_motif.keys()
        sample_pos = sample_profile.keys()

        intersect_pos = set(map(lambda pos: str(math.floor(float(pos))), motif_pos)) & set(map(lambda pos: str(math.floor(float(pos))), sample_pos))

        insert_pattern = r"\d+\.\d+"
        for pos in motif_pos:
            pos_standardised = str(math.floor(float(pos)))
            
            if pos_standardised in Disregard_InDels:
                continue
            
            elif (pos_standardised in intersect_pos) and (re.search(insert_pattern, pos) is None):
                sample_profile[pos] = best_motif[pos]

            elif re.search(insert_pattern, pos) is not None:
                pass




    

