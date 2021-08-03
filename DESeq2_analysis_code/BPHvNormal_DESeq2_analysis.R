# DESeq2 analysis script for prostate BPH v normal prostate
# note: the graphs produced by this script are not in the paper
# this code produces figures 8,9 and 10 in my dissertation


library(BiocManager)
library(DESeq2)
library(readr)
library(dplyr)
library(magrittr)
library(tximport)

# read in normal data run table from SRA
run_table_normal = read_csv('../../normal_bulk/SraRunTable - Copy.txt')

# select useful columns
run_table_normal= select(run_table_normal, Run, BioSample, Type, phenotype)

# read in BPH data run table from SRA
run_table_BPH = read_csv('NCL Project/Progression_data_unprocessed_read_counts/BPHRunTable.txt')

# select useful columns
run_table_BPH= select(run_table_BPH, Run, BioSample, Type, phenotype)

# use dplyr full_join to merge run tables
run_table_BPHvnormal = full_join(run_table_BPH,run_table_normal)


# count files for tximport
BPHvnormal_sample_ids = pull(run_table_BPHvnormal, Run)
BPHvnormal_count_files = paste0('../../BPHCaPRCounts/',BPHvnormal_sample_ids, '/quant.sf')

# set names so the count data is easier to understand
BPHvnormal_conditions = pull(run_table_BPHvnormal, Type)
names(BPHvnormal_count_files) = paste(BPHvnormal_conditions)

# read in mapped transcripts file
transcripts_mapped = read_csv('../Documents/R testing/transcripts_mapped.csv',
                              col_names = c("Transcript_ID","Gene_ID"))

# read in read count data using tximport
BPHvnormal_count_data = tximport(files = BPHvnormal_count_files, type='salmon',
                                ignoreTxVersion = TRUE, tx2gene = transcripts_mapped)


BPHvnormal_deseq_data = DESeqDataSetFromTximport(txi = BPHvnormal_count_data, colData = 
                                                  run_table_BPHvnormal, design = ~Type)

# DESeq analysis
BPHvnormal_deseq_data = DESeq(BPHvnormal_deseq_data)

# get normalized counts
BPHvnormal_counts = counts(BPHvnormal_deseq_data, normalized=TRUE)

# convert matrix to dataframe
BPHvnormal_dataframe = data.frame(BPHvnormal_counts)

# put GeneID into dataframe so it can be searched for
BPHvnormal_dataframe <- tibble::rownames_to_column(BPHvnormal_dataframe, "GeneID")

# plot dispersion estimates
plotDispEsts(BPHvnormal_deseq_data)

# use deseq results function to get results
BPHvnormal_results_table = results(BPHvnormal_deseq_data , contrast = c("Type", "BPH","normal"))

# to data frame to filter for non zeros
BPHvnormal_results_df_gsea = data.frame(counts(BPHvnormal_deseq_data, normalized = TRUE))
BPHvnormal_results_df_full_cases = BPHvnormal_results_df_gsea[!rowSums(BPHvnormal_results_df_gsea[, -1] == 0.0000000) == (ncol(BPHvnormal_results_df_gsea)-1), ]

# write counts to csv for GSEA
write.csv(x=BPHvnormal_results_df_full_cases, file = 'BPHvnormal_deseq_counts.csv', row.names = TRUE)

# check each gene for differential expression
BPHvnormal_results_table$dif_exp = BPHvnormal_results_table$padj < 0.05 & abs(BPHvnormal_results_table$log2FoldChange) > 1
BPHvnormal_dif_exp_results = data.frame(BPHvnormal_results_table)
write.csv(x=BPHvnormal_results_table, file = 'BPHvnormal_total.csv', row.names = TRUE)
# add gene ID to dataframe
BPHvnormal_dif_exp_results$Gene_ID = rownames(BPHvnormal_dif_exp_results)

# use dplyr to filter for differentially expressed genes
BPHvnormal_dif_exp_filtered = filter(BPHvnormal_dif_exp_results, padj < 0.05)
BPHvnormal_dif_exp_filtered_final = filter(BPHvnormal_dif_exp_filtered, abs(log2FoldChange) > 1)

#extract list of differentially expressed genes
BPHvnormal_gene_list_dif_exp = BPHvnormal_dif_exp_filtered_final$Gene_ID
BPHvnormal_gene_list_df = data.frame(BPHvnormal_gene_list_dif_exp)
BPHvnormal_lfc_filtered = BPHvnormal_dif_exp_filtered_final$log2FoldChange
BPHvnormal_gene_lfc = data.frame(BPHvnormal_gene_list_dif_exp, BPHvnormal_lfc_filtered)

# save list to CSV file
row.names(BPHvnormal_gene_lfc) = NULL
write.csv(x=BPHvnormal_gene_lfc, file='BPHvnormal.csv', row.names = FALSE)

# MA plot
plotMA(BPHvnormal_deseq_data , alpha = 0.01, main = "BPH v Normal prostate")

#log transformation
BPHvnormal_rlog_data <- rlogTransformation(BPHvnormal_deseq_data , blind = TRUE)

# PCA plot w/ rlog transformation
plotPCA(BPHvnormal_rlog_data , intgroup = 'Type')
# PCA w/ variance stabilizing
BPHvnormal_var_transform = varianceStabilizingTransformation(BPHvnormal_deseq_data)
plotPCA(BPHvnormal_var_transform , intgroup = 'Type')
plotPCA(BPHvnormal_var_transform , intgroup = 'phenotype')

# heatmap

library('gplots')

BPHvnormal_dist_rl = dist(t(assay(BPHvnormal_rlog_data)))
BPHvnormal_dist_mat = as.matrix(BPHvnormal_dist_rl)
heatmap.2(BPHvnormal_dist_mat , trace = "none", Rowv = FALSE, Colv = FALSE, main = "BPH v Normal prostate")

library("pheatmap")

rld <- normTransform(BPHvnormal_deseq_data)
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
res <- results(BPHvnormal_deseq_data, tidy=TRUE, contrast=c("Type", "BPH","normal")) %>%
  arrange(padj, pvalue) %>%
  tibble::as_tibble()

# write to csv to access row names object easier
write.csv(x = res, file = 'BPHvnormal_resobject.csv', row.names = FALSE)

# select genes of interest
goi <- res$row[c(2322, 13431, 1498, 888)]
goi

# AR
AR = res$row[26446]
tcounts <- t((counts(BPHvnormal_deseq_data[AR, ], normalized=TRUE, replaced=FALSE)+.5)) %>%
  merge(colData(BPHvnormal_deseq_data), ., by="row.names") %>%
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
tcounts <- t(log2((counts(BPHvnormal_deseq_data[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(BPHvnormal_deseq_data), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

ggplot(tcounts, aes(Type, expression, fill = Type,)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  theme_classic() +
  labs(x="Prostate Stage", 
       y="Expression (Log Normalized counts)",
       fill = 'Prostate progression stage',
       title="AR, p53, PTEN and MYC expression (BPH v Normal)")