# Differential Gene Expression Analysis: Tumor vs Normal Tissue Comparison

## Overview
This project performs a comprehensive differential gene expression analysis comparing tumor and normal tissue samples using RNA-Seq data from GEO dataset **GSE50760**. The analysis identifies significantly differentially expressed genes (DEGs), explores their biological functions through enrichment analysis, and investigates protein-protein interaction networks to reveal potential biomarkers and therapeutic targets.

---

## Project Objectives
- Identify genes significantly up- or down-regulated in tumor tissue compared to normal tissue.
- Perform functional enrichment analysis using Gene Ontology (GO), KEGG, and Reactome pathway databases.
- Construct protein-protein interaction (PPI) networks to discover hub genes involved in tumor biology.
- Provide insights into molecular mechanisms underlying tumor development.

---

## Dataset
- **Source:** GEO dataset GSE50760
- **Data Type:** Normalized FPKM RNA-Seq expression values
- **Samples:** Tumor vs Normal tissue samples with phenotype metadata

---

## Methods
- Data preprocessing and quality control (PCA, correlation analysis)
- Differential expression analysis using the **limma** package
- Multiple testing correction with Benjamini-Hochberg FDR
- Functional enrichment analysis via GO, KEGG, and Reactome
- PPI network analysis using STRING database

---

## Key Findings
- Strong immune system dysregulation in tumor samples, particularly adaptive immune responses.
- Disruption of cell cycle regulation and TP53 activity pathways.
- Metabolic reprogramming consistent with cancer hallmarks.
- Identification of potential hub proteins as novel therapeutic targets.

---

## Results Visualization
The repository contains the following key plots to support the analysis:
- Principal Component Analysis (PCA) plot showing sample clustering
- Volcano plot of DEGs highlighting significance and fold change
- Heatmap of top differentially expressed genes
- Barplots for GO and KEGG pathway enrichment
- Protein-protein interaction network visualization

---


