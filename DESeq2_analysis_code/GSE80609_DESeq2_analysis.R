# DESeq2 analysis script for prostate cancer v BPH
# note: the graphs produced by this script are not in the paper
# this script produces figures 11 and 12 in my dissertation

library(BiocManager)
library(DESeq2)
library(readr)
library(dplyr)
library(magrittr)
library(tximport)

# read in run table from SRA
run_table = read_csv('NCL Project/Progression_data_unprocessed_read_counts/SraRunTable.txt')

# select useful columns from run table
run_table_selected = dplyr::select(run_table, Run, BioSample,"Sample Name", 
                            "normal/tumor", progression_step)
run_table_bph = filter(run_table_selected, progression_step == 'BPH')
run_table_CaP = filter(run_table_selected, progression_step == 'CaP')

# use dplyr full_join to merge filtered dfs
run_table_bphvCaP = full_join(run_table_bph,run_table_CaP)

# count files for tximport
bphvCaP_sample_ids = pull(run_table_bphvCaP, Run)
bphvCaP_count_files = paste0('NCL Project/progression_data_processed_rcounts/',bphvCaP_sample_ids, '/quant.sf')

# set names so the count data is easier to understand
bphvCaP_conditions = pull(run_table_bphvCaP, progression_step)
names(bphvCaP_count_files) = paste(bphvCaP_conditions)

# read in mapped transcripts file
transcripts_mapped = read_csv('../Documents/R testing/transcripts_mapped.csv',
                              col_names = c("Transcript_ID","Gene_ID"))

# read in read count data using tximport
bphvCaP_count_data = tximport(files = bphvCaP_count_files, type='salmon',
                              ignoreTxVersion = TRUE, tx2gene = transcripts_mapped)


bphvCaP_deseq_data = DESeqDataSetFromTximport(txi = bphvCaP_count_data, colData = 
                                                run_table_bphvCaP, design = ~progression_step)

# DESeq analysis
bphvCaP_deseq_data = DESeq(bphvCaP_deseq_data)

# get normalized counts
bphvCaP_counts = counts(bphvCaP_deseq_data, normalized=TRUE)

# convert matrix to dataframe
bphvCaP_dataframe = data.frame(bphvCaP_counts)

# put GeneID into dataframe so it can be searched for
bphvCaP_dataframe <- tibble::rownames_to_column(bphvCaP_dataframe, "GeneID")

# plot dispersion estimates
plotDispEsts(bphvCaP_deseq_data)

# use deseq results function to get results
bphvCaP_results_table = results(bphvCaP_deseq_data , contrast = c("progression_step", "CaP","BPH"))

# to data frame to filter for non zeros
bphvCaP_results_df_gsea = data.frame(counts(bphvCaP_deseq_data, normalized = TRUE))
bphvCaP_results_df_full_cases = bphvCaP_results_df_gsea[!rowSums(bphvCaP_results_df_gsea[, -1] == 0.0000000) == (ncol(bphvCaP_results_df_gsea)-1), ]

# write counts to csv for GSEA
write.csv(x=bphvCaP_results_df_full_cases, file = 'CaPvBPH_deseq_counts.csv', row.names = TRUE)

# check each gene for differential expression
bphvCaP_results_table$dif_exp = bphvCaP_results_table$padj < 0.05 & abs(bphvCaP_results_table$log2FoldChange) > 1
bphvCaP_dif_exp_results = data.frame(bphvCaP_results_table)
write.csv(x=bphvCaP_results_table, file = 'CaPvBPH_total.csv', row.names = TRUE)
# add gene ID to dataframe
bphvCaP_dif_exp_results$Gene_ID = rownames(bphvCaP_dif_exp_results)

# use dplyr to filter for differentially expressed genes
bphvCaP_dif_exp_filtered = filter(bphvCaP_dif_exp_results, padj < 0.05)
bphvCaP_dif_exp_filtered_final = filter(bphvCaP_dif_exp_filtered, abs(log2FoldChange) > 1)

#extract list of differentially expressed genes
bphvCaP_gene_list_dif_exp = bphvCaP_dif_exp_filtered_final$Gene_ID
bphvCaP_gene_list_df = data.frame(bphvCaP_gene_list_dif_exp)
bphv_CaP_lfc_filtered = bphvCaP_dif_exp_filtered_final$log2FoldChange
bphvCaP_gene_lfc = data.frame(bphvCaP_gene_list_dif_exp, bphv_CaP_lfc_filtered)

# save list to CSV file
row.names(bphvCaP_gene_lfc) = NULL
write.csv(x=bphvCaP_gene_lfc, file='CaPvBPH.csv', row.names = FALSE)

# MA plot
plotMA(bphvCaP_results_table , alpha = 0.01, main = "GSE80609 Dataset: CaP v BPH")

#log transformation
bphvCaP_rlog_data <- rlogTransformation(bphvCaP_deseq_data , blind = TRUE)

# PCA plot w/ rlog transformation
plotPCA(bphvCaP_rlog_data , intgroup = 'progression_step')
# PCA w/ variance stabilizing
bphvCaP_var_transform = varianceStabilizingTransformation(bphvCaP_deseq_data)
plotPCA(bphvCaP_var_transform , intgroup = 'progression_step')

# heatmap

library('gplots')

bphvCaP_dist_rl = dist(t(assay(bphvCaP_rlog_data)))
bphvCaP_dist_mat = as.matrix(bphvCaP_dist_rl)
heatmap.2(bphvCaP_dist_mat , trace = "none", Rowv = FALSE, Colv = FALSE, main = "GSE80609 Dataset: CaP v BPH")

# annotate genes

library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
bphvCaP_results_table$symbol <- mapIds(org.Hs.eg.db, 
                                       keys=row.names(bphvCaP_results_table), 
                                       column="SYMBOL", 
                                       keytype="ENSEMBL",
                                       multiVals="first")
bphvCaP_results_table$entrez <- mapIds(org.Hs.eg.db, 
                                       keys=row.names(bphvCaP_results_table), 
                                       column="GENENAME", 
                                       keytype="ENSEMBL",
                                       multiVals="first")

library("pheatmap")

rld <- normTransform(bphvCaP_deseq_data)
genes <- c("ENSG00000169083","ENSG00000141510","ENSG00000171862","ENSG00000136997")
mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("progression_step")])
pheatmap(mat, annotation_col = anno, cluster_cols = FALSE, labels_row = c("AR", "p53","PTEN","MYC"),
         main = "GSE80609 Dataset: CaP v BPH")

# boxplot

# load packages
library(tidyr)
library(ggplot2)
library(ggthemes)

# create results object as tibble with tidyr
res <- results(bphvCaP_deseq_data, tidy=TRUE, contrast=c("progression_step", "CaP","BPH")) %>%
  arrange(padj, pvalue) %>%
  tibble::as_tibble()

# write to csv to access row names object easier
write.csv(x = res, file = 'CancervBPH_resobject.csv', row.names = FALSE)

# select genes of interest
goi <- res$row[c(26446, 23161, 21289, 89)]
goi

# AR
AR = res$row[26446]
tcounts <- t((counts(bphvCaP_deseq_data[AR, ], normalized=TRUE, replaced=FALSE)+.5)) %>%
  merge(colData(bphvCaP_deseq_data), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(AR)+1):ncol(.))

ggplot(tcounts, aes(progression_step, expression, fill = progression_step,)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  theme_classic() +
  labs(x="Prostate Stage", 
       y="Expression (Normalized counts)",
       fill = 'Prostate progression stage',
       title="AR expression (Cancer v BPH)")

# goi
tcounts <- t(log2((counts(bphvCaP_deseq_data[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(bphvCaP_deseq_data), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

ggplot(tcounts, aes(progression_step, expression, fill = progression_step,)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  theme_classic() +
  labs(x="Prostate Stage", 
       y="Expression (Log Normalized counts)",
       fill = 'Prostate progression stage',
       title="AR, p53, PTEN and MYC expression (Cancer v BPH)")