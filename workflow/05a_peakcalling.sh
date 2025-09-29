#!/bin/bash
# Step 6: Peak calling with MACS2
# Input: Uniquely mapped, sorted BAM files in ../results/bam/
#        Pattern: <sample>_WT_uniq_sorted.bam and <sample>_inp_uniq_sorted.bam
# Output: Peak files (.narrowPeak, .xls, etc.) in ../results/peaks/
#
# Notes:(Parameters can vary based on the experiment executed)
# - Uses --nomodel and --extsize to improve peak detection
# - Uses --keep-dup all to avoid filtering duplicates
# - Genome size (-g) set to 1.2e7 (approx. S. cerevisiae)

# -----------------------------
# Create output directories
# -----------------------------
mkdir -p ../results/peaks/replicates
mkdir -p ../results/peaks/pulldown

# -----------------------------
# Step 6a. Peak calling per replicate
# -----------------------------
for treat in ../results/bam/*_WT_uniq_sorted.bam
do
    base=$(basename "$treat" _WT_uniq_sorted.bam)
    control=../results/bam/${base}_inp_uniq_sorted.bam

    if [ -f "$control" ]; then
        echo "Calling peaks for $base ..."
        macs2 callpeak \
            -t "$treat" \
            -c "$control" \
            -f BAM -g 1.2e7 \
            --nomodel --extsize 200 --keep-dup all \
            -n ${base}_peaks \
            --outdir ../results/peaks/replicates
    else
        echo "No matching input found for $base, skipping..."
    fi
done

# -----------------------------
# Step 6b. Combined pulldown (all WT vs all Inputs)
# -----------------------------
treatments=$(ls ../results/bam/*_WT_uniq_sorted.bam | tr '\n' ' ')
controls=$(ls ../results/bam/*_inp_uniq_sorted.bam | tr '\n' ' ')

if [ -n "$treatments" ] && [ -n "$controls" ]; then
    echo "Calling peaks on combined replicates (pulldown)..."
    macs2 callpeak \
        -t $treatments \
        -c $controls \
        -f BAM -g 1.2e7 \
        --nomodel --extsize 200 --keep-dup all \
        -n pulldown_all \
        --outdir ../results/peaks/pulldown
else
    echo "Not enough BAM files for pulldown analysis."
fi

echo " Peaks called for individual replicates and pulldown."
