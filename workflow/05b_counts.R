#!/usr/bin/env Rscript
# Quantification using featureCounts
# Input: BAM files from ../results/bam/*_uniq_sorted.bam
# Annotation: GTF file (genes) or BED file (peaks)
# Output: Counts matrix in ../results/counts/
#Note:-will depend upon the gtf annotation file for further downstream analysis
#gene-level counts → use the yeast GTF (S288c_annotation.gtf).
#peak-level counts → convert your MACS2 narrowPeak file(s) into SAF (4-column table) and pass it to annot.ext
# -----------------------------
# Load libraries
# -----------------------------
suppressMessages({
  library(Rsubread)    # featureCounts
})

# -----------------------------
# Create output directory
# -----------------------------
outdir <- "../results/counts"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# -----------------------------
# Input BAM files
# -----------------------------
bam_files <- list.files("../results/bam", pattern = "_uniq_sorted.bam$", full.names = TRUE)

if (length(bam_files) == 0) {
  stop("No BAM files found in ../results/bam/. Run 04_postprocess.sh first.")
}

# -----------------------------
# Annotation file
# -----------------------------
# OPTION 1: Count reads per gene
# gtf_file <- "../reference/S288c_annotation.gtf"

# OPTION 2: Count reads in peaks called by MACS2 (narrowPeak format)
# Convert narrowPeak to SAF (Simplified Annotation Format)
# Example SAF columns: GeneID, Chr, Start, End, Strand

# For now, we assume gene-level counts with a GTF file:
gtf_file <- "../reference/S288c_annotation.gtf"

if (!file.exists(gtf_file)) {
  stop("Annotation file not found: provide a GTF file in ../reference/")
}

# -----------------------------
# Run featureCounts
# -----------------------------
fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "gene",
  GTF.attrType = "ID",
  useMetaFeatures = TRUE,
  nthreads = 4
)

# -----------------------------
# Save outputs
# -----------------------------
# Counts table
counts <- fc$counts
colnames(counts) <- gsub("_uniq_sorted.bam", "", basename(bam_files))

# Write to CSV
write.csv(counts, file.path(outdir, "featureCounts_gene_counts.csv"))

# Log file
sink(file.path(outdir, "featureCounts_summary.txt"))
print(fc$stat)
sink()

cat(" Counts matrix saved in", outdir, "\n")
