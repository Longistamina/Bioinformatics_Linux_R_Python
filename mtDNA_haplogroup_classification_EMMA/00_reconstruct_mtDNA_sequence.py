from Bio import SeqIO
import os
import json
import math
import logging
from argparse import ArgumentParser
import shutil

logger = logging.getLogger(__name__)

class JsonGenerator:
  
    def __init__(self, input_dir, output_dir, ref_path):
        """Initialize JsonGenerator with input and output directories and reference path."""
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.ref_path = ref_path
        
        # Delete output directory if it exists
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        # Create output directory
        os.makedirs(output_dir)

    @staticmethod
    def normalize_pos(pos):
        """Convert position to integer if it's a whole number."""
        if isinstance(pos, int):
            return pos
        elif pos.is_integer():
            return int(pos)
        return pos

    def load_reference_sequence(self):
        """Load reference sequence from FASTA file."""
        for record in SeqIO.parse(self.ref_path, "fasta"):
            return list(record.seq)
        logger.error(f"Failed to load reference sequence from {self.ref_path}")
        return []

    def create_consensus_sequence(self, ref_seq, intervals):
        """Create initial consensus sequence with regions marked."""
        con_seq = list(ref_seq)

        # Mark positions outside of regions as 'N'
        for i in range(len(con_seq)):
            pos = i + 1
            in_region = any(start <= pos <= end for start, end in intervals.values())
            if not in_region:
                con_seq[i] = "N"
                
        return con_seq 
    
    def apply_variants(self, con_seq, variants):
        """Apply SNP and deletion variants to the consensus sequence."""
        for variant in variants:
            pos = variant['pos']
            if variant['ref'] != "-":  # Not an insertion
                idx = int(pos - 1)
                if con_seq[idx] != variant['ref']:
                    logger.warning(f"Different reference base at {pos}, expected {variant['ref']} found {con_seq[idx]}")
                con_seq[idx] = variant['seq']
        
        return con_seq

    def extract_region_sequences(self, ref_seq, con_seq, intervals):
        """Extract sequence regions from reference and consensus sequences."""
        ref_seq_regions = {}
        sample_seq_regions = {}
        
        for region, interval in intervals.items():
            start, end = interval[0], interval[1]
            
            ref_seq_regions[region] = {
                'start': start,
                'end': end,
                'seq': con_seq[start - 1:end]
            }
            
            sample_seq_regions[region] = {
                'start': start,
                'end': end,
                'seq': ref_seq[start - 1:end]
            }
            
        return ref_seq_regions, sample_seq_regions

    def apply_insertions(self, ref_seq, con_seq, variants, ref_seq_regions, sample_seq_regions):
        """Apply insertion variants to sequences."""
        # Sort variants by position
        variants_pos = [self.normalize_pos(variant['pos']) for variant in variants]
        variants_pos.sort()
        
        sorted_variants = {}
        for pos in variants_pos:
            for variant in variants:
                if self.normalize_pos(variant['pos']) == pos:
                    sorted_variants[pos] = variant
                    break
        
        count_insert = 0
        for pos, variant in sorted_variants.items():
            # Process insertion variants
            if variant['ref'] == "-":  # insert condition
                insert_pos = math.floor(float(pos)) + count_insert
                
                # Insert into reference and consensus sequences
                ref_seq.insert(insert_pos, variant['ref'])
                con_seq.insert(insert_pos, variant['seq'])
                
                # Update region-specific sequences
                for region_name in sample_seq_regions.keys():
                    region_start = sample_seq_regions[region_name]['start']
                    region_end = sample_seq_regions[region_name]['end']
                    
                    if region_start <= pos <= region_end:
                        region_insert_pos = math.floor(float(pos)) + count_insert - (region_start - 1)
                        ref_seq_regions[region_name]['seq'].insert(region_insert_pos, variant['ref'])
                        sample_seq_regions[region_name]['seq'].insert(region_insert_pos, variant['seq'])
                
                count_insert += 1
        
        return ref_seq, con_seq, ref_seq_regions, sample_seq_regions

    def gen_seq_from_variants_regions(self, variants, intervals):
        """Generate sequences from variants and regions."""
        # Load reference sequence
        ref_seq = self.load_reference_sequence()
        if not ref_seq:
            return "", "", {}, {}
        
        # Create initial consensus sequence
        con_seq = self.create_consensus_sequence(ref_seq, intervals)
        
        # Apply SNP and deletion variants
        con_seq = self.apply_variants(con_seq, variants)
        
        # Extract region sequences
        ref_seq_regions, sample_seq_regions = self.extract_region_sequences(ref_seq, con_seq, intervals)
        
        # Apply insertion variants
        ref_seq, con_seq, ref_seq_regions, sample_seq_regions = self.apply_insertions(
            ref_seq, con_seq, variants, ref_seq_regions, sample_seq_regions
        )
        
        return "".join(ref_seq), "".join(con_seq), ref_seq_regions, sample_seq_regions

    def write_json_file(self, output_path, save_regions):
        """Write sequences to json format file."""
        with open(output_path, "w") as f:
            output_dict = {}
            for region_name in save_regions.keys():
                
                interval = save_regions[region_name][0] # interval = [start, end]
                seq = save_regions[region_name][1]#.replace("-", "")

                output_dict[region_name] = {
                    "interval": interval,
                    "seq": seq
                }
            # Save output_dir into JSON file:
            json.dump(output_dict, f, indent=4, sort_keys=False)

    def generate_reviewed_json(self):
        """Generate JSON files for all samples in the input directory."""
        success_count = 0
        failure_count = 0
        for json_file in os.listdir(self.input_dir):
            if not json_file.endswith(".json"):
                json_path = os.path.join(self.input_dir, json_file, f'{json_file}.json')
                if not os.path.exists(json_file):
                    logger.warning(f"JSON file not found for sample: {json_file}")
                    failure_count += 1
                    continue
            else:
                json_path = json_path = os.path.join(self.input_dir, json_file)           
            
            try:    
                batch = json.load(open(json_path))

                for sample in batch.keys():
                    
                    sample_data = batch[sample]
                    
                    variants = sample_data['variants']['snps'] + sample_data['variants']['insertions'] + sample_data['variants']['deletions']
                    
                    intervals = sample_data['intervals']
                    for region in intervals.keys():
                        intervals[region] = intervals[region][0] # remove the status of nested list 

                    # Standardize HV2 interval if needed
                    if (intervals['HV2'] != None) and (intervals['HV2'] == [73, 340]):
                        intervals['HV2'] = [68, 340]
                    
                    # Generate sequences
                    insert_ref_seq, delete_con_seq, ref_seq_regions, sample_seq_regions = self.gen_seq_from_variants_regions(variants, intervals)
                    
                    # Prepare sequences for saving
                    save_regions = { 
                        region_name: [intervals[region_name], "".join(sample_seq_regions[region_name]['seq'])] \
                        for region_name in intervals.keys()
                        }
                    # Add deletion and insertion sequences
                    save_regions["delete_in_sample"] = [[1, len(delete_con_seq)], delete_con_seq]
                    save_regions["insert_rCRS"] = [[1, len(insert_ref_seq)], insert_ref_seq]
                    
                    # Write to JSON file
                    output_path = os.path.join(self.output_dir, f"{sample.split('.')[0]}.json")
                    self.write_json_file(output_path, save_regions)
        
                    success_count += 1
                
            except Exception as e:
                logger.error(f"Error processing sample {sample}: {e}")
                failure_count += 1
        
        return success_count, failure_count

def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(description="Generate JSON files from variant data")
    parser.add_argument("-i", "--input", required=True, help="Path to input directory")
    parser.add_argument("-o", "--output", required=True, help="Path to output directory")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA file")
    return parser.parse_args()

def main():
    """Main entry point."""
    args = parse_args()
    
    # Create generator with required reference path
    generator = JsonGenerator(args.input, args.output, args.reference)
    success, failure = generator.generate_reviewed_json()
    
    logger.info(f"JSON generation complete. Successful: {success}, Failed: {failure}")

if __name__ == "__main__":
    main()
