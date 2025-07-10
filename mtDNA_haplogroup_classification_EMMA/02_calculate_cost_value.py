from siuba import *
import numpy as np
import json
import logging
from argparse import ArgumentParser
from multiprocessing import Pool

logger = logging.getLogger(__name__)

class CostValueCalculator:

    def __init__ (self, input_path, output_path, cores = 2):
        '''
        input_path: path to the input json file storing fluctuation rates of each variant position (fluctuation_rate.json)
        output_path: path to the output json file storing cost value of each variant position (cost_value.json)
        cores: number of CPU cores to use for parallel processing
        '''
        self.input_path = input_path
        self.output_path = output_path
        self.cores = cores
    

    def load_fluctuation_rate (self):
        with open(self.input_path, "r") as json_file:
            self.fluctuation_rate = json.load(json_file)
        
        return self.fluctuation_rate
    

    def calculate_cost_value(self, ref_pos):
        dict_cost_value = dict()
        dict_cost_value[ref_pos] = dict()

        for alpha_base in self.fluctuation_rate[ref_pos].keys():
            dict_cost_value[ref_pos][alpha_base] = dict()
            stable = 1 - sum(self.fluctuation_rate[ref_pos][alpha_base].values())
        
            for beta_base, fluctuation in self.fluctuation_rate[ref_pos][alpha_base].items():
                dict_cost_value[ref_pos][alpha_base][beta_base] = np.log10(stable/fluctuation)/3
                    
        return dict_cost_value
    

#-------------------------------------------------------#                
#----------------- main() function ---------------------#
#-------------------------------------------------------#
def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(description="Generate JSON files from variant data")
    parser.add_argument("-i", "--input", required=True, help="Path to input file fluctuation_rate.json")
    parser.add_argument("-o", "--output", required=True, help="Path to output file cost_value.json")
    parser.add_argument("-c", "--cores", type=int, default=2, required=False, help="The number of CPU cores to be used for parallel computing, default is 2. (Example: set -c 8 means using 8 cores to process 8 different samples at a time)")
    return parser.parse_args()

def main():
    # calculator = CostValueCalculator(
    #     input_path = "ref_dir/fluctuation_rate.json", 
    #     output_path = "ref_dir/cost_value.json",  
    #     cores = 14
    # )
    args = parse_args()
    calculator = CostValueCalculator(
        input_path = args.input,
        output_path = args.output,
        cores = args.cores
    )

    try:
        calculator.load_fluctuation_rate()
        
        arg_list = [(ref_pos,) for ref_pos in calculator.fluctuation_rate.keys()]

        with Pool(processes = calculator.cores) as pool:
            result_list = pool.starmap(func = calculator.calculate_cost_value, iterable = arg_list)
        
        cost_value_output_dict = dict()
        for pos_cost_value in result_list:
            cost_value_output_dict.update(pos_cost_value)
        
        with open(calculator.output_path, "w") as json_file:
            json.dump(cost_value_output_dict, json_file, indent = 4)

        logger.info("All positions processed successfully")
    
    except Exception as e:
        logger.error(f"Error in main {e}")


if __name__ == '__main__':
    main()

