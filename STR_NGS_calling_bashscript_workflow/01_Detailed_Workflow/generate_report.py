import numpy as np
import pandas as pd
import glob as glb
import os
import re

np.set_printoptions(suppress = True)
pd.set_option('display.width', 1000)
#pd.set_option('future.no_silent_downcasting', True)

list_df = []

for marker_path in glb.glob("05_calling_annotation/*.txt"):
    
    lines = open(marker_path, 'r').readlines()
    
    n_skip = sum(1 for line in lines if line.startswith('#')) - 1
    calling_result = pd.read_csv(marker_path, sep = '\t', skiprows = n_skip)
    
    marker = os.path.splitext(marker_path.split('/')[1])[0]
    
    marker_type = ''
    for line in lines:
        if line.startswith("#Position"):
            marker_type = "X-STR" if "chrX" in line \
            else "Y-STR" if "chrY" in line \
            else "A-STR"
            break

    calling_result = (
        calling_result
        .rename(columns={"#Sample":"Sample"})
        .assign(Sample=lambda df: df['Sample'].str.replace('_trimmed_aligned_sorted.bam', ''))
        .assign(STR_Type=[marker_type]*calling_result.shape[0])
        .assign(STR_id=[marker]*calling_result.shape[0])
    )

    #Handle heterogeneous allele
    if re.match(".*/.*", calling_result.loc[0, "Allele"]):
        allele_hetero = (
            calling_result
            .pipe(lambda df: df['Allele'].str.extract(r"(.*)/(.*)", expand = True).rename(columns={0:"Allele_1", 1:"Allele_2"}))
            .assign(
                Allele_1=lambda df: pd.to_numeric(df["Allele_1"], errors = 'coerce'),
                Allele_2=lambda df: pd.to_numeric(df["Allele_2"], errors = 'coerce')
                )
            .assign(
                Allele_1=lambda df: df[["Allele_1", "Allele_2"]].min(axis=1).where(df.notna().all(axis=1)),
                Allele_2=lambda df: df[["Allele_1", "Allele_2"]].max(axis=1)
                )
        )

        calling_result_new = (
            pd.concat([calling_result, allele_hetero], axis = 1)
            .pipe(lambda df: df[["Sample", "STR_id", "STR_Type", "Allele_1", "Allele_2", "Probability", "TotalReads", "ValidReads", "SpanReads", "AlleleList"]])
        )
        
    #Handle homogeneous allele
    else:
        allele_homo = (
            calling_result
            .pipe(lambda df: df['Allele'].str.extract(r"(.*\..*)()", expand = True).rename(columns={0:"Allele_1", 1:"Allele_2"}))
            .assign(
                Allele_1=lambda df: pd.to_numeric(df["Allele_1"], errors = 'coerce'),
                Allele_2=pd.Series([np.nan]*calling_result.shape[0])
                )
        )

        calling_result_new = (
            pd.concat([calling_result, allele_homo], axis = 1)
            .pipe(lambda df: df[["Sample", "STR_id", "STR_Type", "Allele_1", "Allele_2", "Probability", "TotalReads", "ValidReads", "SpanReads", "AlleleList"]])
        )

    list_df.append(calling_result_new)

df_final = (
    pd.concat(list_df, axis = 0)
    .sort_values(by=['Sample', 'STR_Type', 'STR_id'])
    )
    
df_final.to_excel("05_calling_annotation/calling_dataframe/str_calling_report.xlsx")
df_final.to_csv("05_calling_annotation/calling_dataframe/str_calling_report.tsv", sep = '\t')

print("\n===>>> STR calling report generated in ./05_calling_annotation/calling_dataframe/ !!!\n")
