# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 17:08:16 2021

@author: o0jl0
"""


# gets the embryonic and cancer common gene list

import pandas as pd

# read in gene sets

normalvSC = pd.read_csv("normalvSC.csv")
BPHvnormal = pd.read_csv("BPHvnormal.csv")
CaPvBPH = pd.read_csv("CaPvBPH.csv")

## first find common genes between all samples

# merge cancer and stem cell dataframes using left join
merged = CaPvBPH.merge(normalvSC, how='inner')

# select gene list column only
merged = merged['gene_list_dif_exp']

# convert series back to dataframe
merged = merged.to_frame()

# merge dataframe again, this time with the BPH dataframe
merged = merged.merge(BPHvnormal, how='inner')

# select gene list column only
merged = merged['gene_list_dif_exp']

# save cancer exclusive genes to csv
merged.to_csv('common_gene_list.csv', index = False)

# merge with original dataframe to get df with log fold change values
common_genes = CaPvBPH.merge(merged, how = 'inner')

common_genes.to_csv('common_genes_ALL.csv', index = False)





## find common genes between cancer and stem cells only

cancer_embryonic_genes = CaPvBPH.merge(normalvSC, how='inner')

# select gene list only
cancer_embryonic_genes = cancer_embryonic_genes['gene_list_dif_exp']

# convert series back to dataframe
cancer_embryonic_genes = cancer_embryonic_genes.to_frame()

# merge two dataframes on left join
cancer_embryonic_genes = cancer_embryonic_genes.merge(common_genes, how='left', indicator=True)

# select rows that are left only
cancer_embryonic_genes = cancer_embryonic_genes[cancer_embryonic_genes['_merge']=='left_only']

# select genes only
cancer_embryonic_genes = cancer_embryonic_genes['gene_list_dif_exp']

# make copy of object for later use
cancer_embryonic_genes_copy = cancer_embryonic_genes

# save gene list to csv
cancer_embryonic_genes.to_csv('cancer_embryonic_genes_gene_list.csv', index = False)

# merge back with original dataframe to get LFC values
cancer_embryonic_genes = CaPvBPH.merge(cancer_embryonic_genes, how = 'inner')

# save object which includes LFC and genes
cancer_embryonic_genes.to_csv('cancer_embryonic_genes_gene_list_LFC.csv', index = False)

# save normalvSC LFC
cancer_embryonic_genes_copy = normalvSC.merge(cancer_embryonic_genes_copy, how = 'inner')

cancer_embryonic_genes_copy.to_csv('normal_LFC_cancer_embryonic_gene_list.csv', index = False)
