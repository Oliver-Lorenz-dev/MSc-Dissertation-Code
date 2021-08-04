# -*- coding: utf-8 -*-
"""
Created on Tue May 18 13:09:10 2021

@author: c0059478
"""

import pandas as pd

'''script extracts gene list of up regulated and downregulated genes \
and writes them into separate files for pathway analysis using DAVID'''


# BPHvCaP

# read in CSV file
CaPvBPH_df = pd.read_csv('CaPvBPH.csv')

# check which genes are upregulated and downregulated
CaPvBPH_df['bphv_CaP_lfc_filtered'] = CaPvBPH_df['bphv_CaP_lfc_filtered'] > 0

# filter for upregulation
CaPvBPH_df_up = CaPvBPH_df[CaPvBPH_df['bphv_CaP_lfc_filtered'] == True]

# filter for downregulation
CaPvBPH_df_down = CaPvBPH_df[CaPvBPH_df['bphv_CaP_lfc_filtered'] == False]

# select gene lists
CaPvBPH_gene_list_up = CaPvBPH_df_up.iloc[:, 0]
CaPvBPH_gene_list_down = CaPvBPH_df_down.iloc[:, 0]

# drop indexes
CaPvBPH_gene_list_up = CaPvBPH_gene_list_up.reset_index(drop = True)
CaPvBPH_gene_list_down = CaPvBPH_gene_list_down.reset_index(drop = True)

# convert series to strings
CaPvBPH_gene_list_up_str = CaPvBPH_gene_list_up.to_string(index = False, header = 'None')
CaPvBPH_gene_list_down_str = CaPvBPH_gene_list_down.to_string(index = False, header = 'None')

# write gene lists to files

with open('gene_list_CaPvBPH_up.txt', 'w') as f:
    f.write(CaPvBPH_gene_list_up_str)
    
with open('gene_list_CaPvBPH_down.txt', 'w') as f:
    f.write(CaPvBPH_gene_list_down_str)

# remove spaces from files

with open('gene_list_CaPvBPH_up.txt', 'r') as f:
    lines = f.readlines()

# remove spaces
lines = [line.replace(' ', '') for line in lines]

with open('gene_list_CaPvBPH_up.txt', 'w') as f:
    f.writelines(lines)
    
with open('gene_list_CaPvBPH_down.txt', 'r') as f:
    lines = f.readlines()

# remove spaces
lines = [line.replace(' ', '') for line in lines]

with open('gene_list_CaPvBPH_down', 'w') as f:
    f.writelines(lines)


# BPHvnormal
    
# read in CSV file
BPHvnormal_df = pd.read_csv('BPHvnormal.csv')

# check which genes are upregulated and downregulated
BPHvnormal_df['BPHvnormal_lfc_filtered'] = BPHvnormal_df['BPHvnormal_lfc_filtered'] > 0

# filter for upregulation
BPHvnormal_df_up = BPHvnormal_df[BPHvnormal_df['BPHvnormal_lfc_filtered'] == True]

# filter for downregulation
BPHvnormal_df_down = BPHvnormal_df[BPHvnormal_df['BPHvnormal_lfc_filtered'] == False]

# select gene lists
BPHvnormal_gene_list_up = BPHvnormal_df_up.iloc[:, 0]
BPHvnormal_gene_list_down = BPHvnormal_df_down.iloc[:, 0]

# drop indexes
BPHvnormal_gene_list_up = BPHvnormal_gene_list_up.reset_index(drop = True)
BPHvnormal_gene_list_down = BPHvnormal_gene_list_down.reset_index(drop = True)

# convert series to strings
BPHvnormal_gene_list_up_str = BPHvnormal_gene_list_up.to_string(index = False, header = 'None')
BPHvnormal_gene_list_down_str = BPHvnormal_gene_list_down.to_string(index = False, header = 'None')

# write gene lists to files

with open('gene_list_BPHvnormal_up.txt', 'w') as f:
    f.write(BPHvnormal_gene_list_up_str)
    
with open('gene_list_BPHvnormal_down.txt', 'w') as f:
    f.write(BPHvnormal_gene_list_down_str)

# remove spaces from files

with open('gene_list_BPHvnormal_up.txt', 'r') as f:
    lines = f.readlines()

# remove spaces
lines = [line.replace(' ', '') for line in lines]

with open('gene_list_BPHvnormal_up.txt', 'w') as f:
    f.writelines(lines)
    
with open('gene_list_BPHvnormal_down.txt', 'r') as f:
    lines = f.readlines()

# remove spaces
lines = [line.replace(' ', '') for line in lines]

with open('gene_list_BPHvnormal_down.txt', 'w') as f:
    f.writelines(lines)    
    
# normalvSC
    
# read in CSV file
normalvSC_df = pd.read_csv('normalvSC.csv')

# check which genes are upregulated and downregulated
normalvSC_df['normalvSC_lfc_filtered'] = normalvSC_df['normalvSC_lfc_filtered'] > 0

# filter for upregulation
normalvSC_df_up = normalvSC_df[normalvSC_df['normalvSC_lfc_filtered'] == True]

# filter for downregulation
normalvSC_df_down = normalvSC_df[normalvSC_df['normalvSC_lfc_filtered'] == False]

# select gene lists
normalvSC_gene_list_up = normalvSC_df_up.iloc[:, 0]
normalvSC_gene_list_down = normalvSC_df_down.iloc[:, 0]

# drop indexes
normalvSC_gene_list_up = normalvSC_gene_list_up.reset_index(drop = True)
normalvSC_gene_list_down = normalvSC_gene_list_down.reset_index(drop = True)

# convert series to strings
normalvSC_gene_list_up_str = normalvSC_gene_list_up.to_string(index = False, header = 'None')
normalvSC_gene_list_down_str = normalvSC_gene_list_down.to_string(index = False, header = 'None')

# write gene lists to files

with open('gene_list_normalvSC_up.txt', 'w') as f:
    f.write(normalvSC_gene_list_up_str)
    
with open('gene_list_normalvSC_down.txt', 'w') as f:
    f.write(normalvSC_gene_list_down_str)
    
# remove spaces from files

with open('gene_list_normalvSC_up.txt', 'r') as f:
    lines = f.readlines()

# remove spaces
lines = [line.replace(' ', '') for line in lines]

with open('gene_list_normalvSC_up.txt', 'w') as f:
    f.writelines(lines)
    
with open('gene_list_normalvSC_down.txt', 'r') as f:
    lines = f.readlines()

# remove spaces
lines = [line.replace(' ', '') for line in lines]

with open('gene_list_normalvSC_down.txt', 'w') as f:
    f.writelines(lines)
    
    
# cancer v normal
    
# read in CSV file
cancervnormal_df = pd.read_csv('cancervnormal.csv')

# check which genes are upregulated and downregulated
cancervnormal_df['cancervnormal_lfc_filtered'] = cancervnormal_df['cancervnormal_lfc_filtered'] > 0

# filter for upregulation
cancervnormal_df_up = cancervnormal_df[cancervnormal_df['cancervnormal_lfc_filtered'] == True]

# filter for downregulation
cancervnormal_df_down = cancervnormal_df[cancervnormal_df['cancervnormal_lfc_filtered'] == False]

# select gene lists
cancervnormal_gene_list_up = cancervnormal_df_up.iloc[:, 0]
cancervnormal_gene_list_down = cancervnormal_df_down.iloc[:, 0]

# drop indexes
cancervnormal_gene_list_up = cancervnormal_gene_list_up.reset_index(drop = True)
cancervnormal_gene_list_down = cancervnormal_gene_list_down.reset_index(drop = True)

# convert series to strings
cancervnormal_gene_list_up_str = cancervnormal_gene_list_up.to_string(index = False, header = 'None')
cancervnormal_gene_list_down_str = cancervnormal_gene_list_down.to_string(index = False, header = 'None')

# write gene lists to files

with open('gene_list_cancervnormal_up.txt', 'w') as f:
    f.write(cancervnormal_gene_list_up_str)
    
with open('gene_list_cancervnormal_down.txt', 'w') as f:
    f.write(cancervnormal_gene_list_down_str)
    
# remove spaces from files

with open('gene_list_cancervnormal_up.txt', 'r') as f:
    lines = f.readlines()

# remove spaces
lines = [line.replace(' ', '') for line in lines]

with open('gene_list_cancervnormal_up.txt', 'w') as f:
    f.writelines(lines)
    
with open('gene_list_cancervnormal_down.txt', 'r') as f:
    lines = f.readlines()

# remove spaces
lines = [line.replace(' ', '') for line in lines]

with open('gene_list_cancervnormal_down.txt', 'w') as f:
    f.writelines(lines)  