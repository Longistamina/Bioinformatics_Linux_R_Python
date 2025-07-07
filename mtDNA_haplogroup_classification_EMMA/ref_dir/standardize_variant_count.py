'''                 WHAT WILL THIS SCRIPT DO?
1. 
'''

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

import pandas as pd
import os

df_raw = pd.read_csv(
    "ref_dir/raw_variant_count.csv", 
    sep = ",", 
    usecols = ["pos", "ref", "alt", "haplogroup", "cnt", "total", "percentage"]
)

df_raw['pos'] = df_raw['pos'].astype(str)

# Find all indices to process
del_indices = df_raw[df_raw['ref'].str.len() > 1].index
insert_indices = df_raw[(df_raw['alt'].str.len() > 1) & (df_raw['alt'] != "del")].index

# Process deletions
single_del_rows = []
for df_index in del_indices:
    for pos_delta, del_base_code in enumerate(df_raw.loc[df_index, "ref"]):
        
        del_pos = str(int(df_raw.loc[df_index, "pos"]) + pos_delta)
        del_base_tuple = IUPAC_CODE.get(del_base_code) # Handle heteroplasmy if exists

        for del_base in del_base_tuple:
            new_single_del_row = pd.DataFrame({
                "pos": [del_pos],
                "ref": [del_base],
                "alt": ["gap"],
                "haplogroup": [df_raw.loc[df_index, "haplogroup"]],
                "cnt": [df_raw.loc[df_index, "cnt"]],
                "total": [df_raw.loc[df_index, "total"]],
                "percentage": [df_raw.loc[df_index, "percentage"]]
            })
            single_del_rows.append(new_single_del_row)

# Process insertions
single_insert_rows = []
for df_index in insert_indices:
    for pos_delta, insert_base_code in enumerate(df_raw.loc[df_index, "alt"]):
        
        if pos_delta == 0:  # Skip reference base
            continue
        
        insert_pos = f'{df_raw.loc[df_index, "pos"]}.{pos_delta}'
        insert_base_tuple = IUPAC_CODE.get(insert_base_code) # Handle heteroplasmy if exists

        for insert_base in insert_base_tuple:
            new_single_insert_row = pd.DataFrame({
                "pos": [insert_pos],
                "ref": ["gap"],
                "alt": [insert_base],
                "haplogroup": [df_raw.loc[df_index, "haplogroup"]],
                "cnt": [df_raw.loc[df_index, "cnt"]],
                "total": [df_raw.loc[df_index, "total"]],
                "percentage": [df_raw.loc[df_index, "percentage"]]
            })
            single_insert_rows.append(new_single_insert_row)

# Combine all indices to remove
all_indices_to_remove = list(set(del_indices) | set(insert_indices))

# Create final processed DataFrame
df_standardized = df_raw.copy().drop(axis = 0, index = all_indices_to_remove)

# Concatenate all parts
all_new_rows = single_del_rows + single_insert_rows
if all_new_rows:  # Only concatenate if there are new rows
    df_standardized = pd.concat([df_standardized] + all_new_rows, axis=0, ignore_index=True)

# Clean up and convert pos to string
df_standardized = (df_standardized
                .drop_duplicates(keep = "first")
                .assign(alt = lambda df: df['alt'].str.replace("del", "gap"))
                .assign(pos = lambda df:  df["pos"].astype(str))  # Quoting the pos column
               )

# Save to CSV
out_csv_path = "ref_dir/standardized_variant_count.csv"
if os.path.exists(out_csv_path):
    os.remove(out_csv_path)
    
df_standardized.to_csv(out_csv_path, index = False)
