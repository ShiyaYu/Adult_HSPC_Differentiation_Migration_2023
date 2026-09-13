# 🩸 Multi-Tissue Hematopoietic Stem and Progenitor Cell (HSPC) Atlas

> **Deciphering the principles of HSPC differentiation, migration, and tissue-specific recruitment at single-cell resolution.**

## 📖 Abstract

Hematopoiesis is a multi-tissue process orchestrated by the coordinated differentiation and migration of hematopoietic stem and progenitor cells (HSPCs). While bone marrow (BM) hematopoiesis is well-characterized, the principles governing HSPC migration across tissues remain unclear. 

To address this, we integrated single-cell transcriptomes of **321,465 human CD34⁺ cells** from bone marrow (BM), peripheral blood (PB), mobilized PB, thymus, and spleen. This unified model of multi-tissue hematopoiesis elucidates the distinct regulatory principles underlying broad systemic dissemination versus niche-specific recruitment.

---

## 📂 Repository Structure

The `codes` directory contains the R scripts used for the bioinformatic analysis, structured sequentially according to the analytical workflow:

```text
📁 codes/
 ├── 📄 01_HSPC_Atlas_Integration.R            # Data loading, quality control, and integration of 321,465 cells across BM, PB, mPB, thymus, and spleen.
 ├── 📄 02_Trajectory_Inference_Monocle3.R     # Pseudotime and developmental trajectory analysis of HSPC differentiation and migration.
 ├── 📄 03_Tissue_specific_DEG_Inference.R     # Identification of differentially expressed genes (DEGs) driving tissue-specific HSPC states.
 ├── 📄 04_Cell_Cell_Communication.R           # Receptor-ligand interaction analysis (e.g., MDK-NCL, CXCL12-CXCR4).
 ├── 📄 05_Thymus_ETP_Analysis.R               # Focused analysis on the ETP_like subset, TSP subsets, and thymic homing mechanisms.
 └── 📄 06_MF_Analysis.R                       # Spatial transcriptomics and functional assay downstream analyses.
