# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 13:11:42 2021

@author: o0jl0
"""

# compare gene sets for venn diagram

import pandas as pd

# read in gene sets

normalvSC = pd.read_csv("normalvSC.csv")
BPHvnormal = pd.read_csv("BPHvnormal.csv")
CaPvBPH = pd.read_csv("CaPvBPH.csv")

# get total sizes of gene lists
print("normal v SC size:", len(normalvSC))
print("BPH v normal size:", len(BPHvnormal)) 
print("Cancer v BPH", len(CaPvBPH))  
# merge dataframes and get value counts

# compare normalvSC and BPHvnormal gene sets 
normalvSC_BPHvnormal = normalvSC.merge(BPHvnormal, how = 'inner')
print('NormalvSC and BPHvNormal')
print(normalvSC_BPHvnormal.value_counts())

# compare normalvSC and CaPvBPH gene sets
normalvSC_CaPvBPH = normalvSC.merge(CaPvBPH, how = 'inner')
print('NormalvSC and CaPvBPH')
print(normalvSC_CaPvBPH.value_counts())

# compare BPHvnormal and CaPvBPH gene sets
BPHvnormal_CaPvBPH = BPHvnormal.merge(CaPvBPH, how = 'inner')
print('BPHvnormal and CaPvBPH')
print(BPHvnormal_CaPvBPH.value_counts())

# compare all conditions to find common genes in all gene sets
normalvSC_BPHvnormal_CaPvBPH = normalvSC_BPHvnormal.merge(CaPvBPH, how = 'inner')
print('iPSC33 + normal + BPH + cancer')
print(normalvSC_BPHvnormal_CaPvBPH.value_counts())

# save merged dataframes
normalvSC_BPHvnormal.to_csv('normalvSC_BPHvnormal.csv',index=False)
normalvSC_CaPvBPH .to_csv('normalvSC_CaPvBPH.csv',index=False)
BPHvnormal_CaPvBPH.to_csv('BPHvnormal_CaPvBPH.csv',index=False)
BPHvnormal_CaPvBPH.to_csv('BPHvnormal_CaPvBPH.csv',index=False)
normalvSC_BPHvnormal_CaPvBPH.to_csv('normalvSC_BPHvnormal_CaPvBPH',index=False)
