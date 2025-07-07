import pandas as pd
from siuba import *


df_variant_count_profile = pd.read_csv(
            "ref_dir/standardized_variant_count.csv", 
            sep = ",", 
            usecols = ["pos", "ref", "alt", "haplogroup", "cnt", "total"]
        ).assign(pos = lambda df: df['pos'].astype(str).str.replace('.0', '', regex = False))


df_alpha_query = (
    df_variant_count_profile
        >> filter(_.pos == "5", _.ref == "A")
        >> mutate(ref_cnt = _.total - _.cnt)
        >> select(~_.alt, ~_.cnt)
        >> rename(alt =_.ref, cnt = _.ref_cnt)
)

print(df_alpha_query)

df_beta_query = (
    df_variant_count_profile
        >> filter(_.pos == "5", _.alt == "C")
        >> select(~_.ref)
)

fill_values = {
    'alt_alpha': "A",
    'cnt_alpha': 0,
    'alt_beta': "C",
    "cnt_beta": 0
}

df_merge = pd.merge(
    left = df_alpha_query, 
    right = df_beta_query, 
    on = ['pos', 'haplogroup', 'total'], 
    how = 'outer', 
    suffixes = ['_alpha', '_beta']
).fillna(value = fill_values)

print(df_merge)