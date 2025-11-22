---
title: "Bulk RNA-seq Normalization Methods"
author: "Nasim Rahmatpour"
date: "2025-11-21"
slug: "bulk-rnaseq-normalization-methods"
categories: ["Bulk RNA-seq"]
tags: ["Normalization", "RNA-seq", "Bioinformatics"]
draft: false
---

Normalization is essential in differential gene expression (DGE) analysis and general data exploration. Raw read counts do not directly reflect true expression levels because they depend on:

• sequencing depth  
• gene length  
• RNA composition of each sample  

Without proper normalization, DGE can produce biased or misleading results.

Below are the most common normalization methods, their intended applications, and whether they are appropriate for DGE.

---

## 1. CPM (Counts Per Million)

• Adjusts only for sequencing depth.  
• Appropriate for comparing the same gene across different samples.  
• Not appropriate for comparing different genes within a sample (does not correct for gene length).  
• ❌ Not suitable for DGE analysis (does not correct compositional bias).

---

## 2. RPKM (Reads Per Kilobase per Million Reads)

• Normalizes for sequencing depth and gene length.  
• Designed for single-end reads.  
• Useful mainly for within-sample comparison of genes.  
• ❌ Not ideal for comparing across samples (library totals differ).  
• ❌ Not suitable for DGE.

---

## 3. FPKM (Fragments Per Kilobase per Million)

• Same as RPKM but for paired-end reads.  
• Corrects for sequencing depth + gene length.  
• Historically popular, now rarely used.  
• ❌ Not appropriate for DGE.

---

## 4. TPM (Transcripts Per Million)

• Also normalizes for gene length and sequencing depth, but the order is reversed compared to RPKM/FPKM.  
• TPM forces the sum of normalized expression to be identical for all samples (1,000,000).  
• This allows you to interpret TPM as the proportion of transcripts mapping to each gene.  
• Good for comparing expression of the same gene across samples.  
• Good for comparing different genes within the same sample.  
• ❌ Not suitable for DGE because it does not correct compositional bias.

---

## Why NONE of these (CPM, TPM, RPKM, FPKM) should be used for DGE

These methods correct for sequencing depth and (sometimes) gene length, but they do **not** correct compositional bias.

Compositional bias occurs when:

• one sample has many highly expressed genes (e.g., hemoglobin, ribosomes),  
  which artificially suppresses the apparent expression of other genes.

Tools such as **DESeq2**, **edgeR**, and **limma/voom** solve this by:

• estimating size factors,  
• adjusting for compositional bias,  
• modeling count data with an appropriate distribution,  
• providing valid DGE results.

(I’ll write a separate post explaining DESeq2 and edgeR.)

---

## Why these methods are not appropriate for single-cell RNA-seq

These normalization methods correct only for sequencing depth and sometimes gene length.  
But single-cell RNA-seq data behave very differently. For example, TPM forces every cell to have the same total expression — which is incorrect for single-cell biology.

### ✔️ In bulk RNA-seq

The total RNA content is averaged across **millions of cells**, so:

• Biological differences (cell size, cycle) wash out.  
• True RNA amount per sample is relatively stable.  
• Extraction efficiency is similar across samples.  
• Most variation comes from sequencing depth.  

So, in bulk RNA-seq, **library size ≈ sequencing depth**.

### ❌ In single-cell RNA-seq

Variation in library size comes from:

**Real biological differences:**  
1. Cell type  
2. Cell size  
3. Cell cycle phase  
4. Activation state  
5. Metabolic activity  

**Technical noise:**  
1. Capture efficiency (5–20%, varies per cell)  
2. Dropouts and zero inflation  
3. mRNA loss during lysis  
4. Sequencing depth differences  

Because of this, the observed library size in a single cell is **not** a reliable measure of true RNA content.

TPM (and similar methods) wrongly force every cell to have the same total expression, which:

• removes true biological differences,  
• amplifies technical noise,  
• distorts downstream analyses (clustering, trajectories, DE).

---

## Summary

Bulk and single-cell data behave very differently:

• **Bulk:** library size ≈ sequencing depth → depth-based normalization works.  
• **Single-cell:** library size = biology + noise + capture efficiency + sequencing depth → depth-based normalization fails.

Thus, CPM/TPM/RPKM/FPKM are useful **abundance metrics**,  
but **NOT valid** for differential expression or single-cell normalization.

---

## Note

TPM is sometimes used in **non-UMI** single-cell protocols such as Smart-seq2/3, where gene length matters.  
However, TPM is **not appropriate** for UMI-based single-cell technologies.

---

🔍 For better visualization examples, see this link:  
*(https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)*

