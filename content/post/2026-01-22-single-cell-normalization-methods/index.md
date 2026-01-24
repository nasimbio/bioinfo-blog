---
title: "Single-Cell Normalization Methods"
author: "Nasim Rahmatpour"
date: "2026-01-22"
slug: "single-cell-normalization-methods"
categories: ["Single-cell RNA-seq"]
tags: ["Normalization", "RNA-seq", "Spatial Transcriptomics", "ATAC-seq", "ADT", "Bioinformatics", "Single-cell"]
draft: false
---

In a previous blog post, I discussed common **bulk RNA-seq normalization** methods. Those approaches (CPM, RPKM/FPKM, TPM) are generally **not appropriate for single-cell RNA-seq**, because in single-cell data the observed library size varies due to a mix of:

- **Technical factors** (sequencing depth, capture efficiency, and many zeros/dropouts)
- **Real biological differences** (cell type, cell size, cell-cycle stage, activation/metabolic state)

In **bulk RNA-seq**, these cell-to-cell differences are averaged across millions of cells, and extraction/capture efficiency tends to be more consistent across samples, so library size is more tightly linked to sequencing depth.

In addition, **TPM includes gene-length correction**, which is largely irrelevant for **UMI-based scRNA-seq**, and TPM-style relative scaling can introduce **compositional bias**.

Below, I’ll walk through several common normalization methods used in **single-cell and spatial workflows** (and other related assays) in **Seurat**, and explain what each method is designed to address.

---

## 1) 🧮 LogNormalize (scRNA/snRNA)

**LogNormalize** is a commonly used normalization method for scRNA/snRNA-seq in Seurat. For each cell, it divides each gene’s raw count by the **total counts (library size)** in that cell, multiplies by a fixed **scale factor** (default = 10,000), and then applies a **log1p transformation**.

### ✅ What problems does it target?
- **Sequencing depth / library size differences:** reduces the effect of cells having different total UMIs, so cells are compared based on **relative abundance** rather than raw counts.
- **Capture efficiency (partially):** because low capture efficiency often results in low total UMIs, rescaling can reduce global differences driven by cell-wide totals, but it doesn’t fix dropouts at the gene level.

### ❌ What does it NOT solve?
- **Dropouts / excess zeros:** it does not model the count-generating process, so it cannot distinguish between a gene that is truly off vs. a gene that is expressed but missed (technical dropout).
- **Mean–variance effects / dominance of highly expressed genes:** even after log normalization, highly expressed genes can still drive much of the variation. This is why workflows typically follow with **HVG selection** and **scaling/centering** before PCA.

---

## 2) 📊️ SCTransform (model-based normalization + variance stabilization)

**SCTransform** is a **model-based normalization** and **variance-stabilization** method (commonly used for scRNA and also widely used in spatial transcriptomics workflows). Instead of simply dividing by total counts, it fits a **regularized Negative Binomial (NB) regression** for each gene across cells. The expected count depends on the cell’s depth (typically **total UMI**) and the gene’s baseline behavior learned across cells.

From this fitted model, it computes **Pearson residuals**, which behave roughly like **z-scores**: they reflect how much the observed count deviates from what is expected **after accounting for depth**, scaled by expected variability.

**In simpler terms:** SCTransform asks:  
> Given this cell’s depth and this gene’s typical behavior, what count should I expect?
> Then it converts the observed value into a *standardized deviation* from that expectation.

### ✅ What problems does it target?
- **Sequencing depth / global library size effects:** depth is modeled directly, so counts are interpreted relative to what is expected for that cell’s total UMI.
- **Capture efficiency (partially)::** if a cell has low capture efficiency, it often has low total UMIs; modeling depth reduces the impact of that global shift.
- **Mean–variance effects:** residuals are variance-stabilized, so highly expressed genes don’t automatically dominate PCA just because they have larger raw variance.
- **Many zeros / dropouts (partially):** SCTransform does not remove zeros, but it makes embeddings less sensitive to random dropouts:
  - a zero that is **expected** (low gene + low depth) is treated as **not informative**
  - a zero that is **unexpected** (high depth + usually detected) is treated as **more surprising**

---

### ⚠️ Note (important tradeoff)
These normalization methods involve a tradeoff: they reduce technical variation driven by library size, but they can also **attenuate real biology** correlated with total UMIs (e.g., cell size, activation, or cell-cycle state).

SCTransform additionally uses a model to downweight “expected zeros,” but it still cannot perfectly label a zero as **technical dropout** vs. **true biological absence** for a single gene in a single cell.

---

Next, let’s briefly cover two approaches commonly used for **non-RNA assays** in Seurat/Signac-style workflows:

- **TF-IDF** for scATAC-seq  
- **CLR** for ADT (CITE-seq)

---

## 3) 🧬 TF-IDF (scATAC-seq)

**TF-IDF** is applied to **peak × cell** accessibility matrices. These matrices are typically much sparser than scRNA-seq data and often behave close to **binary** (a peak is detected or not detected in a given cell).

### What it does
- **TF (term frequency):** normalizes each cell by its total number of accessible fragments/peaks → controls **per-cell depth (library size)**.
- **IDF (inverse document frequency):** downweights peaks accessible in many cells (ubiquitous “open everywhere” peaks) and upweights peaks accessible in fewer cells (more informative peaks).
- The TF-IDF matrix is then used as input for **SVD/LSI**, producing a low-dimensional embedding for clustering and visualization.

### ✅ What problems does it target?
- **Per-cell depth / library size differences**
- **Dominance of ubiquitous peaks** that provide little clustering signal
- Makes downstream embedding more robust to sparsity by emphasizing informative non-zero signals (without imputation)

### ❌ What does it NOT solve?
- TF-IDF does **not** impute missing peaks and cannot determine whether a zero is due to technical dropout or truly closed chromatin. It reweights the observed signals to improve downstream embedding.

---

## 4) 🧪 CLR (Centered Log Ratio) for ADT (CITE-seq)

**CLR** is commonly used for **ADT** (protein) counts. ADT behaves differently from RNA: it often has large dynamic range, background signal, and strong **compositional effects** (if one protein is very high, others can look lower simply because totals are constrained).

### How it works (conceptually)
- Add a **pseudocount** (e.g., +1) so zeros don’t break the log transform.
- For each cell, compute a baseline across proteins using the **geometric mean**.
- Convert each protein value into a **log-ratio** relative to that baseline.

So CLR turns ADT into **relative enrichment** values rather than raw counts.

### ✅ What problems does it target?
- **Compositional bias (main target)**
- **Global per-cell scaling / depth effects (partially)**
- Helps stabilize dynamic range for downstream PCA/WNN

### ❌ What does it NOT solve?
- **Ambient/background antibody signal:** CLR does not explicitly subtract background (methods like **DSB** are better for that).
- CLR does not model dropout like SCTransform does; a zero remains a zero (after pseudocount it becomes small).

---

## 📚 Further Reading

- This paper https://pmc.ncbi.nlm.nih.gov/articles/PMC8238499/#sec3 by **Hao et al. (2021)**,  includes practical workflow descriptions for preprocessing and normalization across modalities (e.g., log-normalization or SCTransform for RNA, CLR for ADT, and TF-IDF + LSI for ATAC).

---


