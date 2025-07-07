from siuba import *
import numpy as np
import pandas as pd
import json
import os
import logging
from argparse import ArgumentParser
from multiprocessing import Pool

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

logger = logging.getLogger(__name__)

VARIANT_SET = {"A", "C", "G", "T", "gap"}

class FluctuationRateCalculator:

    def __init__ (self, input_path, output_path, cores = 2):
        '''
        input_path: path to the input csv file storing variants' count by haplogroup (standardized_variant_count.csv)
        output_path: path to the output json file storing fluctuation rates of each variant position (fluctuation_rate.json)
        cores: number of CPU cores to use for parallel processing
        '''
        self.input_path = input_path
        self.output_path = output_path
        self.cores = cores
    

    def load_variant_count_profile (self):
        self.df_variant_count_profile = pd.read_csv(
            self.input_path, 
            sep = ",", 
            usecols = ["pos", "ref", "alt", "haplogroup", "cnt", "total"]
        ).assign(pos = lambda df: df['pos'].astype(str).str.replace('.0', '', regex = False))

        return self.df_variant_count_profile


    def calculate_fluctuation_rate(self, ref_pos):
        fluctuation_rate = dict()
        fluctuation_rate[ref_pos] = dict()

        ref_base = self.df_variant_count_profile[self.df_variant_count_profile['pos'] == ref_pos]['ref'].unique()[0]

        alpha_list = list(VARIANT_SET)       

        for alpha_base in alpha_list:
            fluctuation_rate[ref_pos][alpha_base] = dict()
            
            if alpha_base == ref_base:
                df_alpha_query = (
                    self.df_variant_count_profile
                        >> filter(_.pos == ref_pos, _.ref == ref_base)
                        >> mutate(ref_cnt = _.total - _.cnt)
                        >> select(~_.alt, ~_.cnt)
                        >> rename(alt =_.ref, cnt = _.ref_cnt)
                )
            else:
                df_alpha_query = (
                    self.df_variant_count_profile
                        >> filter(_.pos == ref_pos, _.alt == alpha_base)
                        >> select(~_.ref)
                )

            beta_list = list(VARIANT_SET - {alpha_base})
            
            for beta_base in beta_list:
                if beta_base == ref_base:
                    df_beta_query = (
                        self.df_variant_count_profile
                            >> filter(_.pos == ref_pos, _.ref == ref_base)
                            >> mutate(ref_cnt = _.total - _.cnt)
                            >> select(~_.alt, ~_.cnt)
                            >> rename(alt =_.ref, cnt = _.ref_cnt)
                    )
                else:
                    df_beta_query = (
                        self.df_variant_count_profile
                            >> filter(_.pos == ref_pos, _.alt == beta_base)
                            >> select(~_.ref)
                    )

                fill_values = {
                    'alt_alpha': alpha_base,
                    'cnt_alpha': 0,
                    'alt_beta': beta_base,
                    "cnt_beta": 0
                }

                df_merge = pd.merge(
                    left = df_alpha_query, 
                    right = df_beta_query, 
                    on = ['pos', 'haplogroup', 'total'], 
                    how = 'outer', 
                    suffixes = ['_alpha', '_beta']
                ).fillna(value = fill_values)

                beta_fr = np.minimum(df_merge['cnt_alpha'], df_merge['cnt_beta']).sum() / df_merge['total'].sum()
                
                if (beta_fr == 0) or np.isnan(beta_fr):
                    transition_check = ({alpha_base, beta_base} == {"A", "G"}) or ({alpha_base, beta_base} == {"C", "T"})                  
                    match transition_check:
                        case True:
                            beta_fr = 10**(-6)
                        case _:
                            beta_fr = 10**(-9)
                
                fluctuation_rate[ref_pos][alpha_base][beta_base] = beta_fr
        
        return fluctuation_rate


#-------------------------------------------------------#                
#----------------- main() function ---------------------#
#-------------------------------------------------------#
def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(description="Generate JSON files from variant data")
    parser.add_argument("-i", "--input", required=True, help="Path to input file standardized_variant_count.csv")
    parser.add_argument("-o", "--output", required=True, help="Path to output file fluctuation_rate.json")
    parser.add_argument("-c", "--cores", required = False, help="The number of CPU cores to be used for parallel computing, default is 2. (Example: set -c 8 means using 8 cores to process 8 different samples at a time)")
    return parser.parse_args()

def main():
    # calculator = FluctuationRateCalculator(
    #     input_path = "ref_dir/standardized_variant_count.csv", 
    #     output_path = "ref_dir/fluctuation_rate.json",  
    #     cores = 14
    # )

    args = parse_args()
    calculator = FluctuationRateCalculator(
        input_path = args.input,
        output_path = args.output,
        cores = args.cores
    )

    try:
        calculator.load_variant_count_profile()

        arg_list = [(ref_pos,) for ref_pos in calculator.df_variant_count_profile['pos'].unique()]

        with Pool(processes = calculator.cores) as pool:
            result_list = pool.starmap(func = calculator.calculate_fluctuation_rate, iterable = arg_list)
        
        fluctuation_rate_output_dict = dict()
        for pos_fluctuation_rate in result_list:
            fluctuation_rate_output_dict.update(pos_fluctuation_rate)
        
        with open(calculator.output_path, "w") as json_file:
            json.dump(fluctuation_rate_output_dict, json_file, indent = 4)

        logger.info("All positions processed successfully")
    
    except Exception as e:
        logger.error(f"Error in main {e}")


if __name__ == '__main__':
    main()