# 🩸 Multi-Tissue Hematopoietic Stem and Progenitor Cell (HSPC) Atlas

> **Deciphering the Cross-Tissue Migration Landscape of Hematopoietic Stem and Progenitor Cells at Single-Cell Resolution**

## 📖 Abstract

Hematopoiesis is a multi-tissue process orchestrated by the coordinated differentiation and migration of hematopoietic stem and progenitor cells (HSPCs). While bone marrow (BM) hematopoiesis is well-characterized, the principles governing HSPC migration across tissues remain unclear. To address this, we integrated single-cell transcriptomes of 321,465 human CD34⁺ cells from BM, peripheral blood (PB), mobilized PB, thymus, and spleen. This enabled the construction of a comprehensive, multi-tissue HSPC atlas. We discovered that most HSPC subsets in BM have transcriptional counterparts in PB, revealing a broad migration spectrum from BM to circulation. We further identified key signaling axes (e.g., MDK-NCL, CXCL12-CXCR4) regulating HSPC retention versus mobilization, where high expression of MDK enhanced HSPCs proliferation and inhibited HSPCs migration. Spatial and functional analyses revealed that HSPCs with low CXCR4 expression reside proximal to sinusoids and exhibit higher migratory potential, suggesting a priming state for egress. Focusing on thymic homing, we defined a HSPC subset with expression profile similar to the earliest thymic progenitors (ETP), thus was called ETP_like. ETP_like was clustered into thymus seeding progenitor (TSP) subsets and ETP subsets. We found the concentration gradient of chemokines centered around the medulla plays an important role in lymphoid precursor cells homing to thymus. Analyses of spatial transcriptomics data and chemotaxis assay showed that signals including CCL19-CCR7 are critical for directing TSP recruitment to the thymus. Collectively, this study establishes a unified model of multi-tissue hematopoiesis and elucidates the distinct regulatory principles underlying broad systemic dissemination versus niche-specific recruitment.

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
