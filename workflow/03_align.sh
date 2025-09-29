#!/bin/bash
# Step 3: Alignment of trimmed reads using Bowtie
# Input: Trimmed FASTQ files in ../results/trimmed/
# Reference: S288c genome (S288C_reference_genome_R64-1-1_20110203.tgz),
#            with headers renamed to chr1..chr16, chrmt
# Output: SAM files in ../results/aligned/

# -----------------------------
# Create output directory
# -----------------------------
mkdir -p ../results/aligned
mkdir -p ../reference

# -----------------------------
# Step 3a. Index the reference genome (only if not already done)
# -----------------------------
if [ ! -f ../reference/WT_index.1.bt2 ]; then
    echo "Indexing reference genome..."
    bowtie2-build ../reference/S288c_reference_genome.fasta ../reference/WT_index
else
    echo "Reference index already exists. Skipping indexing."
fi

# -----------------------------
# Step 3b. Align each trimmed FASTQ
# -----------------------------
for file in ../results/trimmed/*_trimmed.fastq.gz
do
    base=$(basename "$file" _trimmed.fastq.gz)
    echo "Aligning $file ..."
    
    bowtie2 --no-unal -p 4 \
        -x ../reference/WT_index \
        -U "$file" \
        --rg-id "$base" --rg "SM:$base" \
        -S ../results/aligned/${base}.sam
done

echo "SAM files are saved in ../results/aligned/"

