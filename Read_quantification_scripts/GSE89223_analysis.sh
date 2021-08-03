#!/usr/bin/bash

# Tumour versus normal pre processing and quantification script

# need to activate conda environment in terminal

# make directories

mkdir QC

mkdir fastpout

cd fastpout

mkdir QC

cd ..

# run fastqc

fastqc *fastq.gz -o QC

# run multiqc

multiqc QC

# run fastp over dataset

for index in {783..814}
do
       fastp -i SRR4453"${index}".fastq.gz -o fastpout/SRR4453"${index}".fastq.gz -R "SRR4453"${index}"" -h "SRR4453"${index}".html"
done

# run fastqc on processed data

cd fastpout

fastqc *fastq.gz -o QC

# run multiqc

multiqc QC

cd ..

# run salmon over fastp processed dataset


for index in {783..814}
do
        salmon quant -i human_v37_index -l A -r fastpout/SRR4453"${index}".fastq.gz --validateMappings -o fastpout/SRR4453"${index}"
done