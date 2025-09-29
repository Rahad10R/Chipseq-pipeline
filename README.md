# 🧬 ChIP-Seq Analysis Pipeline 

![Conda](https://img.shields.io/badge/Conda-ready-green?logo=anaconda)
![R](https://img.shields.io/badge/R-4.2%2B-blue?logo=r)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)

A **reproducible Next-Generation Sequencing (NGS) pipeline** for **ChIP-seq data analysis**.  
Developed as part of my thesis work on **Msh5 protein binding in *Saccharomyces cerevisiae*** (WT vs mutant).  
This pipeline demonstrates **end-to-end analysis**: from raw FASTQ files to **peaks, counts, and plots** — highlighting my skills in **bioinformatics, NGS workflows, and reproducible research**.

---

## 📌 Features

- 🔎 **Quality Control**: FastQC for raw reads  
- ✂️ **Read Trimming**: Adapter & quality trimming with Trimmomatic  
- 🎯 **Alignment**: Bowtie2 mapping to yeast S288c genome  
- 🧹 **Post-processing**: SAM→BAM conversion, filtering unique reads, sorting & indexing (samtools)  
- 📈 **Peak Calling**: MACS2 for individual replicates & combined analysis  
- 📊 **Read Counting**: featureCounts for genes/peaks  
- 🎨 **Visualization**: R scripts for PCA, correlation heatmap, and peak distribution  

---



