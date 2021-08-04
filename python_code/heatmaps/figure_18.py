# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 18:10:55 2021

@author: o0jl0
"""


# heatmap of top embryonic cancer genes

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv('Final_Cancer_EMT_genes.csv', index_col=0)


# make heatmap
fig = plt.figure(figsize=(2,7))
ax = sns.heatmap(data, cmap='RdYlGn_r')
    
# set title
ax.set_title('Expression levels of key embryonic genes in PCa v BPH', fontweight = 'bold')