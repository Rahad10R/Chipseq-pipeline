#!/usr/bin/env Rscript
# Some preliminiary analysis of ChIP-seq data
# Baesd on the count matrix and peakcalling data..
# Note:- Further downstream analysis can be developed based on aim of the experiments
# Inputs:
#   - ../results/counts/featureCounts_gene_counts.csv (from 05_counts.R)
#   - ../results/peaks/replicates/*.narrowPeak (from 06_peakcalling.sh)
#
# Outputs:
#   - Plots in ../results/plots/

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# -----------------------------
# Create output directory
# -----------------------------
plot_dir <- "../results/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# -----------------------------
# Load counts
# -----------------------------
counts_file <- "../results/counts/featureCounts_gene_counts.csv"
if (!file.exists(counts_file)) {
  stop("Counts file not found: Run 05_counts.R first.")
}

counts <- read.csv(counts_file, row.names = 1)
cat("Loaded counts matrix with", nrow(counts), "features and", ncol(counts), "samples.\n")

# -----------------------------
# DESeq2 normalization
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = data.frame(condition = rep("WT", ncol(counts))),
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# -----------------------------
# PCA plot
# -----------------------------
rld <- rlog(dds)
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
ggsave(file.path(plot_dir, "PCA_plot.png"), p, width = 6, height = 5)

# -----------------------------
# Heatmap of sample correlations
# -----------------------------
sample_cor <- cor(norm_counts)
pheatmap(sample_cor, 
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         main = "Sample correlation heatmap",
         filename = file.path(plot_dir, "sample_correlation_heatmap.png"))

# -----------------------------
# Peak counts per replicate (from MACS2 narrowPeak files)
# -----------------------------
peak_dir <- "../results/peaks/replicates"
if (dir.exists(peak_dir)) {
  peak_files <- list.files(peak_dir, pattern = "*.narrowPeak", full.names = TRUE)
  if (length(peak_files) > 0) {
    peak_counts <- sapply(peak_files, function(f) nrow(read.table(f)))
    peak_df <- data.frame(Sample = names(peak_counts), Peaks = peak_counts)

    p2 <- ggplot(peak_df, aes(x = Sample, y = Peaks, fill = Sample)) +
      geom_col() +
      theme_minimal() +
      labs(title = "Number of peaks per replicate", y = "Peaks", x = "Replicate")
    ggsave(file.path(plot_dir, "peaks_per_replicate.png"), p2, width = 6, height = 5)
  } else {
    message("No narrowPeak files found for peak count barplot.")
  }
}

cat("Demo Plots saved in", plot_dir, "\n")
