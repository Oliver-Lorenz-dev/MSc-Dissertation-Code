#!/usr/bin/bash

# prostate cancer v BPH pre processing and quantification script

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

for index in {43..87}
do
       fastp -i SRR34370"${index}"_1.fastq.gz -I SRR34370"${index}"_2.fastq.gz -o fastpout/SRR34370"${index}"_1.fastq.gz -O fastpout/SRR34370"${index}"_2.fastq.gz --detect_adapter_for_pe -c -R "SRR34370"${index}"" -h "SRR34370"${index}".html"
       rm SRR34370"${index}"_1.fastq.gz
       rm SRR34370"${index}"_2.fastq.gz
done

# run fastqc on processed data

cd fastpout

fastqc *fastq.gz -o QC

# run multiqc

multiqc QC

cd ..

# run salmon over fastp processed dataset


for index in {43..87}
do
        salmon quant -i human_v37_index -l A -1 fastpout/SRR34370"${index}"_1.fastq.gz -2 fastpout/SRR34370"${index}"_2.fastq.gz --validateMappings -o fastpout/SRR34370"${index}"
done