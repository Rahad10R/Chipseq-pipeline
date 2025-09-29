#!/bin/bash
# Quality Control for raw FASTQ files using FastQC

# Input: FASTQ files in example_data/
# Output: QC reports in results/qc/

mkdir -p ../results/qc

for file in ../example_data/*.fastq.gz
do
    echo "Running FastQC on $file"
    fastqc $file -o ../results/qc/
done

