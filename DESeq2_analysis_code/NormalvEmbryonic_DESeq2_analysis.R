# DESeq2 analysis script for prostate normal prostate v embryonic (iPSC33)
# note: the graphs produced by this script are not in the paper
# this script produces figures 5,6 and 7 in my dissertation

library(BiocManager)
library(DESeq2)
library(readr)
library(dplyr)
library(magrittr)
library(tximport)

# read in run table from SRA
run_table = read_csv('../../normal_bulk/SraRunTable.txt')

# select useful columns from run table
run_table_normalvSC= select(run_table, Run, BioSample, Type, phenotype, source_name)

# count files for tximport
normalvSC_sample_ids = pull(run_table_normalvSC, Run)
normalvSC_count_files = paste0('../../normal_bulk/',normalvSC_sample_ids, '/quant.sf')

# set names so the count data is easier to understand
normalvSC_conditions = pull(run_table_normalvSC, Type)
names(normalvSC_count_files) = paste(normalvSC_conditions)

# read in mapped transcripts file
transcripts_mapped = read_csv('../Documents/R testing/transcripts_mapped.csv',
                              col_names = c("Transcript_ID","Gene_ID"))

# read in read count data using tximport
normalvSC_count_data = tximport(files = normalvSC_count_files, type='salmon',
                             ignoreTxVersion = TRUE, tx2gene = transcripts_mapped)


normalvSC_deseq_data = DESeqDataSetFromTximport(txi = normalvSC_count_data, colData = 
                                               run_table_normalvSC, design = ~Type)

# DESeq analysis
normalvSC_deseq_data = DESeq(normalvSC_deseq_data)

# get normalized counts
normalvSC_counts = counts(normalvSC_deseq_data, normalized=TRUE)

# convert matrix to dataframe
normalvSC_dataframe = data.frame(normalvSC_counts)

# put GeneID into dataframe so it can be searched for
normalvSC_dataframe <- tibble::rownames_to_column(normalvSC_dataframe, "GeneID")

# plot dispersion estimates
plotDispEsts(normalvSC_deseq_data)

# use deseq results function to get results
normalvSC_results_table = results(normalvSC_deseq_data , contrast = c("Type", "normal","Stem_cell"))

# to data frame to filter for non zeros
normalvSC_results_df_gsea = data.frame(counts(normalvSC_deseq_data, normalized = TRUE))
normalvSC_results_df_full_cases = normalvSC_results_df_gsea[!rowSums(normalvSC_results_df_gsea[, -1] == 0.0000000) == (ncol(normalvSC_results_df_gsea)-1), ]

# write counts to csv for GSEA
write.csv(x=normalvSC_results_df_full_cases, file = 'normalvSC_deseq_counts.csv', row.names = TRUE)

# check each gene for differential expression
normalvSC_results_table$dif_exp = normalvSC_results_table$padj < 0.05 & abs(normalvSC_results_table$log2FoldChange) > 1
normalvSC_dif_exp_results = data.frame(normalvSC_results_table)
write.csv(x=normalvSC_results_table, file = 'normalvSC_total.csv', row.names = TRUE)
# add gene ID to dataframe
normalvSC_dif_exp_results$Gene_ID = rownames(normalvSC_dif_exp_results)

# use dplyr to filter for differentially expressed genes
normalvSC_dif_exp_filtered = filter(normalvSC_dif_exp_results, padj < 0.05)
normalvSC_dif_exp_filtered_final = filter(normalvSC_dif_exp_filtered, abs(log2FoldChange) > 1)

#extract list of differentially expressed genes
normalvSC_gene_list_dif_exp = normalvSC_dif_exp_filtered_final$Gene_ID
normalvSC_gene_list_df = data.frame(normalvSC_gene_list_dif_exp)
normalvSC_lfc_filtered = normalvSC_dif_exp_filtered_final$log2FoldChange
normalvSC_gene_lfc = data.frame(normalvSC_gene_list_dif_exp, normalvSC_lfc_filtered)

# save list to CSV file
row.names(normalvSC_gene_lfc) = NULL
write.csv(x=normalvSC_gene_lfc, file='normalvSC.csv', row.names = FALSE)

# MA plot
plotMA(normalvSC_deseq_data , alpha = 0.01, main = "Normal prostate v iPSC33")

#log transformation
normalvSC_rlog_data <- rlogTransformation(normalvSC_deseq_data , blind = TRUE)

# PCA plot w/ rlog transformation
plotPCA(normalvSC_rlog_data , intgroup = 'Type')
# PCA w/ variance stabilizing
normalvSC_var_transform = varianceStabilizingTransformation(normalvSC_deseq_data)
plotPCA(normalvSC_var_transform , intgroup = 'Type')
plotPCA(normalvSC_var_transform , intgroup = 'phenotype')

# heatmap

library('gplots')

normalvSC_dist_rl = dist(t(assay(normalvSC_rlog_data)))
normalvSC_dist_mat = as.matrix(normalvSC_dist_rl)
heatmap.2(normalvSC_dist_mat , trace = "none", Rowv = FALSE, Colv = FALSE, main = "Pro iPSC33 v GSE80609 Dataset: BPH samples")

library("pheatmap")

rld <- normTransform(normalvSC_deseq_data)
genes <- c("ENSG00000169083","ENSG00000141510","ENSG00000171862","ENSG00000136997")
mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("Type")])
pheatmap(mat, annotation_col = anno, cluster_cols = FALSE, labels_row = c("AR", "p53","PTEN","MYC"))

# boxplot

# load packages
library(tidyr)
library(ggplot2)
library(ggthemes)

# create results object as tibble with tidyr
res <- results(normalvSC_deseq_data, tidy=TRUE, contrast=c("Type", "normal","Stem_cell")) %>%
  arrange(padj, pvalue) %>%
  tibble::as_tibble()

# write to csv to access row names object easier
write.csv(x = res, file = 'NormalvSC_resobject.csv', row.names = FALSE)

# select genes of interest
goi <- res$row[c(15852, 10130, 2327, 2991)]
goi

# AR
AR = res$row[26446]
tcounts <- t((counts(normalvSC_deseq_data[AR, ], normalized=TRUE, replaced=FALSE)+.5)) %>%
  merge(colData(normalvSC_deseq_data), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(AR)+1):ncol(.))

ggplot(tcounts, aes(Type, expression, fill = Type,)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  theme_classic() +
  labs(x="Prostate Stage", 
       y="Expression (Normalized counts)",
       fill = 'Prostate progression stage',
       title="AR expression (Cancer v BPH)")

# goi
tcounts <- t(log2((counts(normalvSC_deseq_data[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(normalvSC_deseq_data), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

ggplot(tcounts, aes(Type, expression, fill = Type,)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  theme_classic() +
  labs(x="Prostate Stage", 
       y="Expression (Log Normalized counts)",
       fill = 'Prostate progression stage',
       title="AR, p53, PTEN and MYC expression (Normal v iPSC33)")