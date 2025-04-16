## 20.440 Project – Lauren Hoorens van Heyningen & Laura Williams

Project Title: Analyzing subsets of dendritic cells involved in the development of peripheral tolerance in the gut using single cell transcriptomics

## Overview and Purpose

This repository contains code and data used to generate Figure 5 from our project report draft, which shows the transcriptomic features of dendritic cells (DCs) involved in peripheral regulatory T cell (pTreg) induction, with a focus on RORγt⁺ antigen-presenting cells (APCs) identified from a single-cell RNA-seq dataset published by Rodrigues et al. (2025, Cell).
More specifically, Figure 5 presents the results of a Gene Ontology (GO) enrichment analysis on cluster 6 — a transcriptionally distinct population of RORγt⁺ DCs — to identify biological processes that may underlie their role in T cell differentiation and oral tolerance induction. The repository contains the scRNAseq data in a separate file that can be downloaded to
the user's personal computer, along with the R script to generate the plots as we did in R Studio. Lastly, this README.md file can be found in the repository as you just did! Whoo!

The overall purpose of this repository is to provide the user to download the relevant raw scRNAseq data, the R script we used, and instructions on how to replicate Figure 5 of our report.

## Data Description:
- Raw count matrices were processed using the Seurat pipeline in R.
- Quality control, clustering, and marker gene detection were performed to identify the RORγt⁺ DC cluster (cluster 6).
- Genes upregulated in cluster 6 were subjected to GO enrichment analysis.

## Folder Structure
We created a folder that contains three important files of the scRNAseq data (barcodes, features, and matrix) in the folder (`scRNAseq_Data_Files/`). Here you can find:
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`
  
## Figure Replication Instructions

In order to replicate the figure, follow these steps:
1. Download the three files from scRNAseq_Data_Files onto personal computer into a folder named 20.440-Project located on your computer's Desktop.
   You can choose to name your folder differently or choose a different location to store these files, but remember to update line 8 of the "Figure 4-5.R" file.
2. Download the "Figure 4-5.R" file from the repository as well, and open it in R Studio.
3. Run the R script in R studio, ensuring that it can access the three scRNAseq data files.
4. The plot from Figure 5 should appear in the "Plots" tab on the bottom right. It may take some time for it to plot!

## Citations for core datasets and methods:
- Rodrigues et al. (2025), *Cell*: [https://doi.org/10.1016/j.cell.2025.03.020](https://doi.org/10.1016/j.cell.2025.03.020)

## Sources:
Single-cell RNA sequencing data from wild-type mice was downloaded from GEO:  
Accession GSE289268 (Rodrigues et al., 2025)

