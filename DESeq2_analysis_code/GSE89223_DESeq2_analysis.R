# DESeq2 analysis script for prostate cancer v normal prostate
# note: the graphs produced by this script are not in the paper
# this script produces figures 14,15 and 16 in my dissertation

library(BiocManager)
library(DESeq2)
library(readr)
library(dplyr)
library(magrittr)
library(tximport)

# read in run table from SRA
run_table = read_csv('../../Russian data/SraRunTable.txt')

# select useful columns from run table
run_table_selected = dplyr::select(run_table, Run, BioSample, Gleason_sum, 
                                   patient_identifier, Tissue)
run_table_cancer = filter(run_table_selected, Tissue == 'tumorous prostate tissue')
run_table_normal = filter(run_table_selected, Tissue == 'normal prostate tissue')

# use dplyr full_join to merge filtered dfs
run_table_cancervnormal = full_join(run_table_cancer,run_table_normal)

# count files for tximport
cancervnormal_sample_ids = pull(run_table_cancervnormal, Run)
cancervnormal_count_files = paste0('../../Russian data/',cancervnormal_sample_ids, '/quant.sf')

# set names so the count data is easier to understand
cancervnormal_conditions = pull(run_table_cancervnormal, Tissue)
names(cancervnormal_count_files) = paste(cancervnormal_conditions)

# read in mapped transcripts file
transcripts_mapped = read_csv('../Documents/R testing/transcripts_mapped.csv',
                              col_names = c("Transcript_ID","Gene_ID"))

# read in read count data using tximport
cancervnormal_count_data = tximport(files = cancervnormal_count_files, type='salmon',
                              ignoreTxVersion = TRUE, tx2gene = transcripts_mapped)


cancervnormal_deseq_data = DESeqDataSetFromTximport(txi = cancervnormal_count_data, colData = 
                                                run_table_cancervnormal, design = ~Tissue)

# DESeq analysis
cancervnormal_deseq_data = DESeq(cancervnormal_deseq_data)

# get normalized counts
cancervnormal_counts = counts(cancervnormal_deseq_data, normalized=TRUE)

# convert matrix to dataframe
cancervnormal_dataframe = data.frame(cancervnormal_counts)

# put GeneID into dataframe so it can be searched for
cancervnormal_dataframe <- tibble::rownames_to_column(cancervnormal_dataframe, "GeneID")

# plot dispersion estimates
plotDispEsts(cancervnormal_deseq_data)

# use deseq results function to get results
cancervnormal_results_table = results(cancervnormal_deseq_data , contrast = c("Tissue", "tumorous prostate tissue","normal prostate tissue"))

# to data frame to filter for non zeros
cancervnormal_results_df_gsea = data.frame(counts(cancervnormal_deseq_data, normalized = TRUE))
cancervnormal_results_df_full_cases = cancervnormal_results_df_gsea[!rowSums(cancervnormal_results_df_gsea[, -1] == 0.0000000) == (ncol(cancervnormal_results_df_gsea)-1), ]

# write counts to csv for GSEA
write.csv(x=cancervnormal_results_df_full_cases, file = 'cancervnormal_deseq_counts.csv', row.names = TRUE)

# check each gene for differential expression
cancervnormal_results_table$dif_exp = cancervnormal_results_table$padj < 0.05 & abs(cancervnormal_results_table$log2FoldChange) > 1
cancervnormal_dif_exp_results = data.frame(cancervnormal_results_table)
write.csv(x=cancervnormal_results_table, file = 'cancervnormal_total.csv', row.names = TRUE)
# add gene ID to dataframe
cancervnormal_dif_exp_results$Gene_ID = rownames(cancervnormal_dif_exp_results)

# use dplyr to filter for differentially expressed genes
cancervnormal_dif_exp_filtered = filter(cancervnormal_dif_exp_results, padj < 0.05)
cancervnormal_dif_exp_filtered_final = filter(cancervnormal_dif_exp_filtered, abs(log2FoldChange) > 1)

#extract list of differentially expressed genes
cancervnormal_gene_list_dif_exp = cancervnormal_dif_exp_filtered_final$Gene_ID
cancervnormal_gene_list_df = data.frame(cancervnormal_gene_list_dif_exp)
cancervnormal_lfc_filtered = cancervnormal_dif_exp_filtered_final$log2FoldChange
cancervnormal_gene_lfc = data.frame(cancervnormal_gene_list_dif_exp, cancervnormal_lfc_filtered)

# save list to CSV file
row.names(cancervnormal_gene_lfc) = NULL
write.csv(x=cancervnormal_gene_lfc, file='cancervnormal.csv', row.names = FALSE)

# MA plot
plotMA(cancervnormal_results_table , alpha = 0.01, main = "GSE89223 Dataset: Cancer v Normal")

#log transformation
cancervnormal_rlog_data <- rlogTransformation(cancervnormal_deseq_data , blind = TRUE)

# PCA plot w/ rlog transformation
plotPCA(cancervnormal_rlog_data , intgroup = 'Tissue')
# PCA w/ variance stabilizing
cancervnormal_var_transform = varianceStabilizingTransformation(cancervnormal_deseq_data)
plotPCA(cancervnormal_var_transform , intgroup = 'Tissue')

# heatmap

library('gplots')

cancervnormal_dist_rl = dist(t(assay(cancervnormal_rlog_data)))
cancervnormal_dist_mat = as.matrix(cancervnormal_dist_rl)
heatmap.2(cancervnormal_dist_mat , trace = "none", Rowv = FALSE, Colv = FALSE, main = "GSE89223 Dataset: Cancer v Normal")