from Bio import SeqIO
import numpy as np
import json
import os
import logging
from argparse import ArgumentParser
import shutil
import re
from multiprocessing import Pool
from termcolor import colored, cprint

logging.basicConfig(level=logging.INFO)
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

class InputStandardizer():
    def __init__(self, input_dir, output_dir, ref_path, cores = 2):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.ref_path = ref_path
        self.cores = cores
        self.dict_sample_profile = dict()

        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir)

    def standardize_existed_input(self, json_file):
        # self.dict_sample_profile = {}       
        # for json_file in os.listdir(self.input_dir):
        #     if not json_file.endswith(".json"):
        #         json_path = os.path.join(self.input_dir, json_file, f'{json_file}.json')
        #         if not os.path.exists(json_path):
        #             logger.warning(f"JSON file not found for sample: {json_file}")
        #             continue
        #     else:
        #         json_path = os.path.join(self.input_dir, json_file) 

        #     try:
        #         json_dict = json.load(open(json_path, "r"))

        json_path = os.path.join(self.input_dir, json_file)
        json_dict = json.load(open(json_path, "r"))
            
        for sample in json_dict.keys():
            self.dict_sample_profile[sample] = dict()

            ranges = ""
            for interval in json_dict[sample]["intervals"].values():
                interval = interval[0]
                ranges += f"{interval[0]}-{interval[1]} "
            self.dict_sample_profile[sample]["ranges"] = ranges.rstrip(' ')

            for type in json_dict[sample]['variants'].keys():
            
                for variant_info in json_dict[sample]['variants'][type]:
                    pos = str(variant_info['pos'])
                    ref = variant_info['ref']
                    alt = variant_info['seq']

                    self.dict_sample_profile[sample][pos] = dict()
            
                    match type:
                        case "deletions":
                            self.dict_sample_profile[sample][pos].update({"ref": ref, "alt": "gap"})
                        case "insertions":
                            self.dict_sample_profile[sample][pos].update({"ref": "gap", "alt": alt})
                        case _:
                            self.dict_sample_profile[sample][pos].update({"ref": ref, "alt": alt})
        
        return self.dict_sample_profile


    def create_new_input(self):
        for record in SeqIO.parse(self.ref_path, "fasta"):
            rCRS = list(record.seq)

        cprint("\nFailed to read input file.", "black", "on_red", ["bold"])
        
        start_message = colored(
            text = "\nWould you like to create new inputs manually? [y/n]: ",
            color = "black",
            on_color = "on_blue",
            attrs = ["bold"]
        )
        yes_no_error = colored(
            text = '>>>InputError: The response must be either "y" or "n", please try again!',
            color = "black",
            on_color = "on_red",
            attrs = ["bold"]
        )

        sample_message = colored(
            text = "\nYour sample name: ",
            color = "black",
            on_color = "on_green",
            attrs = ["bold"]
        )
        sample_error = colored(
            text = '>>>InputError: Your sample name must not be left blank, please try again!',
            color = "black",
            on_color = "on_red",
            attrs = ["bold"]
        )
        
        pos_message = colored(
            text = "\nPosition of mutation: ",
            color = "black",
            on_color = "on_light_yellow",
            attrs = ["bold"]
        )
        pos_type_error = colored(
            text = '>>>InputError: Your pos must be either integer or real numbers, please try again!',
            color = "black",
            on_color = "on_red",
            attrs = ["bold"]
        )
        pos_range_error = colored(
            text = '>>>InputError: Your pos must lie between 1 and 16569, please try again!',
            color = "black",
            on_color = "on_red",
            attrs = ["bold"]
        )

        mutation_message = colored(
            text = "\nMutation: ",
            color = "black",
            on_color = "on_light_cyan",
            attrs = ["bold"]
        )
        mutation_iupac_error = colored(
            text = f">>>InputError: Your mutation must be one of the IUPAC codes {[dna_code for dna_code in IUPAC_CODE.keys()]}, please try again!",
            color = "black",
            on_color = "on_red",
            attrs = ["bold"]
        )

        while True:
            try:
                start_program = input(start_message).lower()
                assert (start_program == "y") or (start_program == "n")
            
            except AssertionError:
                logger.error(yes_no_error)
                continue
            
            else:
                if start_program == "n":
                    break
                else:
                    logger.info(f"\nNOTE 1: your input mutation position must lie between 1 and 16569.\n")
                    logger.info(f"\nNOTE 2: your input mutation must be one of the IUPAC codes {[dna_code for dna_code in IUPAC_CODE.keys()]}.\n")
                    for dna_code, meaning in IUPAC_CODE.items():
                        print(f"{dna_code}: {meaning}")
                    
                    ### Handle sample name input ##
                    while True:
                        try:
                            sample_name = input(sample_message)
                            assert (sample_message != "") and (re.search(r"^\s+\b", sample_message) is None)
                        except AssertionError:
                            logger.error(sample_error)
                            continue

                        else:
                            self.dict_sample_profile[sample_name] = dict()
                            
                            ### Handle position input ##
                            while True:
                                try:
                                    pos = eval(input(pos_message))
                                    assert (pos >= 1) and (pos <= 16569)
                                except Exception:
                                    logger.error(pos_type_error)
                                    continue
                                except AssertionError:
                                    logger.error(pos_range_error)
                                    continue
                                else:
                                    pos = str(pos).replace(".0", "")
                                    self.dict_sample_profile[sample_name][pos] = dict()
                                    
                                    if re.search(r"\d+\.{1}\d+", pos) is not None:
                                        ref = "gap"
                                    else:
                                        ref = rCRS[int(pos) - 1]
                                    self.dict_sample_profile[sample_name][pos]["ref"] = ref

                                    ### Handle mutation input ##
                                    while True:
                                        try:
                                            mutation = input(mutation_message)

                                            mutation_overlap_error = colored(
                                                text = f'>>>InputError: Your mutation overlaps the rCRS reference base "{ref}", please try again!',
                                                color = "black",
                                                on_color = "on_red",
                                                attrs = ["bold"]
                                            )

                                            assert mutation in IUPAC_CODE.keys(), mutation_iupac_error
                                            assert mutation != ref, mutation_overlap_error
                                        
                                        except AssertionError as error_message:
                                            logger.error(error_message)
                                            continue

                                        else:
                                            self.dict_sample_profile[sample_name][pos]["alt"] = mutation
                                            break
                                
                                ### Ask for adding new mutation in the same sample ##
                                next_mutation_message = colored(
                                    text = f'\nWould you like to add new mutation for this "{sample_name}" sample? [y/n]:',
                                    color = "black",
                                    on_color = "on_light_magenta",
                                    attrs = ["bold"]
                                )
                                while True:
                                    try:
                                        next_mutation = input(next_mutation_message).lower()
                                        assert (next_mutation == "y") or (next_mutation == "n")
                                    except AssertionError:
                                        logger.error(yes_no_error)
                                        continue
                                    else:
                                        break
                                if next_mutation == "y":
                                    continue
                                else:
                                    break
                            

                            ### Handle mutation ranges for easier cross-check ##
                            if len(self.dict_sample_profile[sample_name]) == 1:
                                pos = list(self.dict_sample_profile[sample_name].keys())[0]
                                ranges = f"{np.floor(float(pos))-1}-{np.floor(float(pos))+1}"
                            else:
                                pos_sorted = sorted([np.floor(float(pos)) for pos in self.dict_sample_profile[sample_name].keys()])
                                ranges = f"{pos_sorted[0]}-{pos_sorted[1]}"

                            self.dict_sample_profile[sample_name]["ranges"] = ranges.replace(".0", "")
                        
                        
                        ### Ask for adding new sample ##
                        next_sample_message = colored(
                                text = f'\nWould you like to add new sample? [y/n]:',
                                color = "black",
                                on_color = "on_light_magenta",
                                attrs = ["bold"]
                            )
                        while True:
                            try:
                                next_sample = input(next_sample_message).lower()
                                assert (next_sample == "y") or (next_sample == "n")
                            except AssertionError:
                                logger.error(yes_no_error)
                                continue
                            else:
                                break

                        if next_sample == "y":
                            continue
                        else:
                            out_message = colored(
                                    text = "\nInput process ended!",
                                    color = "black",
                                    on_color = "on_green",
                                    attrs = ["bold"]
                                )
                            logger.info(out_message)
                            break
            break
        
        return self.dict_sample_profile


#----------------------------------------#
#----------- main function --------------#
#----------------------------------------#

def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(description="Generate standardised JSON files from variant data")
    parser.add_argument("-i", "--input", required=True, help="Path to input directory")
    parser.add_argument("-o", "--output", required=True, help="Path to output directory")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA file")
    parser.add_argument("-c", "--cores", type=int, default=2, required=False, help="The number of CPU cores to be used for parallel computing, default is 2. (Example: set -c 8 means using 8 cores to process 8 different samples at a time)")
    return parser.parse_args()

def main():
    # input_creator = InputStandardizer(
    # input_dir = "input_dir/",
    # output_dir = "input_dir/standardized/",
    # ref_path = "ref_dir/rcrs.fasta",
    # cores = 10
    # )
    args = parse_args()

    input_creator = InputStandardizer(
    input_dir = args.input,
    output_dir = args.output,
    ref_path = args.reference,
    cores = 10
    )

    try:
        json_file_list = [(json_file,) for json_file in os.listdir(input_creator.input_dir) if json_file.endswith(".json")]
        
        with Pool(processes = input_creator.cores) as pool:
            result_list = pool.starmap(func = input_creator.standardize_existed_input, iterable = json_file_list)
        
        for standardized_sample in result_list:
            input_creator.dict_sample_profile.update(standardized_sample)
    
    except Exception as e:
        error_message = colored(text = f"\n{e}", color = "red", attrs = ["bold", "underline"])
        logger.error(error_message)
        input_creator.create_new_input()
    
    if len(input_creator.dict_sample_profile) > 0:
        for sample, mutation_info in input_creator.dict_sample_profile.items():
            out_json_path = os.path.join(input_creator.output_dir, f"{sample}.json")
            with open(out_json_path, "w") as json_file:
                json.dump(mutation_info, json_file, indent = 4)
        
        end_message = colored(
            text = f"\nAll standardized input files are saved in {input_creator.output_dir}",
            color = "black",
            on_color = "on_green",
            attrs = ["bold"]
        )
        logger.info(end_message)
    
    else:
        end_message = colored(
            text = f"\nNot sample files were created!!!",
            color = "black",
            on_color = "on_yellow",
            attrs = ["bold"]
        )
        logger.warning(end_message)


if __name__ == "__main__":
    main()


# python3 03_standardize_sample_profile.py -i input_dir/ -o input_dir/standardized -r ref_dir/rcrs.fasta -c 10