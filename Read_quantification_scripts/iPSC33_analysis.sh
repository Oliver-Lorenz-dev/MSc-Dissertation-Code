#!/usr/bin/bash

# pre processing and quantification script for iPSC data from the Heer lab

# ACTIVATE CONDA ENV IN TERMINAL

# make directories

mkdir QC

mkdir fastpout

cd fastpout

mkdir QC

cd ..

# fastqc

fastqc *fastq.gz -o QC

# run fastp

fastp -i index254_TCCTGAGC-AGAGTAGA_L001-L002_R1_001.fastq.gz -I index254_TCCTGAGC-AGAGTAGA_L001-L002_R2_001.fastq.gz -o fastpout/index254_TCCTGAGC-AGAGTAGA_L001-L002_R1_001.fastq.gz -O fastpout/index254_TCCTGAGC-AGAGTAGA_L001-L002_R2_001.fastq.gz --detect_adapter_for_pe -c -R "iPSC fastp report" -h "iPSC_fastp_report.html"

# fastqc

cd fastpout

fastqc *fastq.gz -o QC

cd ..

# run salmon

salmon quant -i human_v37_index -l A -1 index254_TCCTGAGC-AGAGTAGA_L001-L002_R1_001.fastq.gz -2 index254_TCCTGAGC-AGAGTAGA_L001-L002_R2_001.fastq.gz --validateMappings -o salmon_iPSC

salmon quant -i human_v37_index -l A -1 fastpout/index254_TCCTGAGC-AGAGTAGA_L001-L002_R1_001.fastq.gz -2 fastpout/index254_TCCTGAGC-AGAGTAGA_L001-L002_R2_001.fastq.gz --validateMappings -o fastpout/salmon_iPSC