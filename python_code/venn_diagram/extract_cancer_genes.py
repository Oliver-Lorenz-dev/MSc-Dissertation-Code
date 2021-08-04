# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 19:29:01 2021

@author: o0jl0
"""

# gets the cancer exclusive gene list

import pandas as pd

# read in gene sets

normalvSC = pd.read_csv("normalvSC.csv")
BPHvnormal = pd.read_csv("BPHvnormal.csv")
CaPvBPH = pd.read_csv("CaPvBPH.csv")

# merge cancer and stem cell dataframes using left join
merged = CaPvBPH.merge(normalvSC, how='left', indicator=True)

# select left rows that are left only
merged = merged[merged['_merge']=='left_only']

# select gene list column only
merged = merged['gene_list_dif_exp']

# convert series back to dataframe
merged = merged.to_frame()

# merge dataframe again, this time with the BPH dataframe
merged = merged.merge(BPHvnormal, how='left', indicator=True)

# select rows that are left only
merged = merged[merged['_merge']=='left_only']

# select gene list column only
merged = merged['gene_list_dif_exp']

# save cancer exclusive genes to csv
merged.to_csv('cancer_exclusive_gene_list.csv', index = False)

# merge with original dataframe to get df with log fold change values
final_cancer_gene_df = CaPvBPH.merge(merged, how = 'inner')

final_cancer_gene_df.to_csv('cancer_genes_lfc.csv', index=False)
