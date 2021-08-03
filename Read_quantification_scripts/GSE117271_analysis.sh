#!/usr/bin/bash

# normal prostate data pre processing and quantification script

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

for index in {606..621}
do
       fastp -i SRR7533"${index}".fastq.gz -o fastpout/SRR7533"${index}".fastq.gz -R "SRR7533"${index}"" -h "SRR7533"${index}".html"
done

# run fastqc on processed data

cd fastpout

fastqc *fastq.gz -o QC

# run multiqc

multiqc QC

cd ..

# run salmon over fastp processed dataset


for index in {606..621}
do
        salmon quant -i human_v37_index -l A -r fastpout/SRR7533"${index}".fastq.gz --validateMappings -o fastpout/SRR7533"${index}"
done