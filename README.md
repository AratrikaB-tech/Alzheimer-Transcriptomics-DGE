Transcriptomic Analysis of Alzheimer's Disease (GSE53697)

---Project Overview---
This project performs a Differential Gene Expression (DGE) analysis on human brain tissue samples (Post-mortem Cerebral Cortex) to identify molecular biomarkers associated with Alzheimer's Disease. Using a cohort of 17 samples (9 Alzheimer's vs. 8 Healthy Controls), I identified a 484-gene signature that characterizes the neurodegenerative state.

---Key Findings---
 * Statistical Power: Identified 484 genes with a raw p-value < 0.05.
 * Clustering: Unsupervised hierarchical clustering (Heatmap) successfully separated AD patients from the Control group based on gene expression profiles.
 * Biological Context: Significant genes are linked to neuroinflammation and metabolic shifts in the brain, relevant to concepts of cellular signaling.
   
---Results Visualization---
Tools, Packages, and Softwares Used
This analysis was conducted using the R Programming Language (R Console) and the following bioinformatics ecosystem:
| R (CRAN) | Environment | Base R console for script execution and data processing. |
| DESeq2 | Bioconductor | Primary engine for differential expression using shrinkage estimation. |
| pheatmap | CRAN | Generation of hierarchical clustered heatmaps. |
| org.Hs.eg.db | Annotation | Mapping Ensembl/Symbol IDs to human genome nomenclature. |
| ggplot2 | Visualization | Custom plotting for Volcano and MA plots. |

---Methodology---
 * Data Acquisition: Raw count extraction from NCBI GEO (GSE53697).
 * Preprocessing: Metadata alignment and Variance Stabilizing Transformation (VST).
 * Statistical Testing: Applied the Wald test for significance via DESeq2.
 * Visualization: Generated Heatmaps for sample-to-sample comparison using Euclidean distance clustering.
