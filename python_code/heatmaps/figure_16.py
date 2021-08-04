# -*- coding: utf-8 -*-
"""
Created on Wed May 12 15:26:52 2021

@author: c0059478
"""

# heatmap of key developmental genes in PCa v normal prostate

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# function which creates heatmap of gene expression of key developmental genes

def check_gene_expression(df):
    
    # read in expression datafile
    df = pd.read_csv(df, index_col=0)
    
    # select lfc column
    df = df['log2FoldChange']
    
    # select genes of interest
    genes_of_interest = df.loc[['ENSG00000142208','ENSG00000169083',
                                'ENSG00000125378','ENSG00000101144',
                                'ENSG00000039068','ENSG00000070193',
                                'ENSG00000129514','ENSG00000125798',
                                'ENSG00000111087','ENSG00000074047',
                                'ENSG00000159184','ENSG00000106211',
                                'ENSG00000136997',
                                'ENSG00000167034','ENSG00000141510',
                                'ENSG00000073282','ENSG00000171862',
                                'ENSG00000146374','ENSG00000164690',
                                'ENSG00000105329','ENSG00000184012',
                                'ENSG00000112715','ENSG00000085741',
                                'ENSG00000002745','ENSG00000114251',
                                'ENSG00000136244','ENSG00000253293',
                                'ENSG00000106031','ENSG00000128714',
                                'ENSG00000075891','ENSG00000107779',
                                'ENSG00000183691','ENSG00000105989',
                                'ENSG00000147655','ENSG00000168036',
                                'ENSG00000205213','ENSG00000104332',
                                'ENSG00000156076','ENSG00000121879',
                                'ENSG00000125398','ENSG00000090776',
                                'ENSG00000115594','ENSG00000115758',
                                'ENSG00000105894','ENSG00000087245',
                                'ENSG00000173083','ENSG00000137573']]
    
    # convert series back to dataframe
    genes_of_interest = genes_of_interest.to_frame()
    
    # map genes to useful identifier
    gene_list = ['AKT','AR','BMP4','BMP7',
                 'E-Cadherin','FGF10','FOXA1','FOXA2','GLI1','GLI2',
                 'HOXB13','HSP27','MYC',
                 'NKX3.1','p53','p63','PTEN','R-spondin3','SHH','TGFB',
                 'TMPRSS2','VEGF','WNT11','WNT16','WNT5a', 'IL-6','HOXA10',
                 'HOXA13','HOXD13','PAX2','BMPR1A','NOG','WNT2','R-spondin2',
                 'CTNNB1','LGR4','SFRP1','WIF1','PI3K','SOX9','EFNB1','IL1R1',
                 'ODC1','PTN','MMP2','HSPE','SULF1']
    
    # add normal gene names to dataframe
    genes_of_interest['gene_name'] = gene_list
    
    # drop ensembl gene id
    genes_of_interest.reset_index(drop=True, inplace=True)
    
    # rename gene_name column for graph
    genes_of_interest.rename(columns = {'gene_name' : 'Gene'}, inplace=True)
    
    # set gene_name column as index
    genes_of_interest.set_index('Gene', inplace=True)
    
    # sort lfc in descending order
    genes_of_interest.sort_values(by = 'log2FoldChange', ascending=False, inplace=True)
    
    # make heatmap
    fig = plt.figure(figsize=(2,10))
    ax = sns.heatmap(genes_of_interest, cmap='RdYlGn_r')
    
    # have to change title manually
    ax.set_title('Expression levels of key genes in PCa v Normal Prostate', fontweight = 'bold')
    
check_gene_expression('cancervnormal_total.csv')