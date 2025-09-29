#!/bin/bash
# Adapter trimming using Trimmomatic

# Input: FASTQ files in example_data/
# Output: Trimmed FASTQ files in results/trimmed/

#Paramaters adjusted for SE pairs, with adaptors file .fa based on the QC analysis of fq files. Also with different read length(QC_report/sample)
# -- vary the MINLEN parameter as well
mkdir -p ../results/trimmed

for file in ../example_data/*.fastq.gz
do
    base=$(basename $file .fastq.gz)
    echo "Trimming $file"
    trimmomatic SE -threads 4 \
        $file ../results/trimmed/${base}_trimmed.fastq.gz \
        ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
