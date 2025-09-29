#!/bin/bash
# Step 4: Post-processing of SAM files using SAMtools
# Author: Rahad Rehas
# Date: Sept 2025
#
# Input: SAM files in ../results/aligned/
# Output: BAM files (raw, uniq, sorted, indexed) in ../results/bam/

# -----------------------------
# Create output directory
# -----------------------------
mkdir -p ../results/bam

# -----------------------------
# Step 4a. Convert SAM → BAM
# -----------------------------
for sam in ../results/aligned/*.sam
do
    base=$(basename "$sam" .sam)
    echo "Converting $sam to BAM..."
    samtools view -bS "$sam" > ../results/bam/${base}.bam
done

# -----------------------------
# Step 4b. Coverage statistics
# Just to Check the Coverage, also adivised to check after trimming to make sure large amounts of reads havent been trimmed off causing a bias in the downstream analysis  -----------------------------
echo "Sample_Name,Total_Reads,Mapped_Reads" > ../results/bam/mapping_stats.csv

for bam in ../results/bam/*.bam
do
    base=$(basename "$bam" .bam)

    # Total reads (all, including unmapped)
    total=$(samtools view -c "$bam")

    # Mapped (primary aligned) reads only (-F 260 removes unmapped + secondary alignments)
    mapped=$(samtools view -c -F 260 "$bam")

    echo "$base,$total,$mapped" >> ../results/bam/mapping_stats.csv
done

echo "Mapping stats written to ../results/bam/mapping_stats.csv"

# -----------------------------
# Step 4c. Extract uniquely mapped reads
# -----------------------------
for bam in ../results/bam/*.bam
do
    base=$(basename "$bam" .bam)
    echo "Extracting uniquely mapped reads for $bam ..."
    samtools view -h -F 4 "$bam" | grep -v "XS:" | samtools view -bS - > ../results/bam/${base}_uniq.bam
done

# -----------------------------
# Step 4d. Sort and index BAM files
# -----------------------------
for bam in ../results/bam/*_uniq.bam
do
    base=$(basename "$bam" .bam)
    echo "Sorting and indexing $bam ..."
    samtools sort "$bam" -o ../results/bam/${base}_sorted.bam
    samtools index ../results/bam/${base}_sorted.bam
done

echo "Step 4 complete: processed BAMs are ready in ../results/bam/"
