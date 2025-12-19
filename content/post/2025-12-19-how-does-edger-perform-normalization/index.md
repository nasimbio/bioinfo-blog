---
title: "How Does edgeR Perform Normalization?"
author: "Nasim Rahmatpour"
date: 2025-12-19
slug: "edgeR-normalization-explained"
draft: false
categories: ["RNA-seq", "Normalization"]
tags: ["Normalization", "RNA-seq", "DESeq2", "DGE", "RNA-seq", "bulk-RNA"]
math: true
---

## How Does edgeR Perform Normalization?

In my previous blog post, I explained why **library size** and **compositional bias** must be accounted for in differential gene expression (DGE) analysis.  
**edgeR**, similar to **DESeq2**, is a widely used tool for DGE analysis. In this post, I explain **how edgeR performs normalization step by step**, focusing on its default method: **TMM (Trimmed Mean of M values)**.

---

## Step-by-step normalization in edgeR

### 1. Remove genes with zero expression in all samples

Genes that have zero counts across all samples do not contain useful information and would lead to undefined log ratios. These genes are therefore excluded from the normalization procedure.

---

### 2. Pick one sample as the reference

A reference sample is selected and used to normalize all other samples against it.

A **good reference sample** is the most *average* sample, meaning its expression distribution is representative of the dataset.  
A **bad reference sample** would be an extreme one, such as a sample with very small library size or one dominated by a few highly expressed genes.

In practice, edgeR selects a reference sample through the following steps:

- Scale each sample by its total library size (divide each gene count by the total counts of that sample).
- For each sample, calculate the **75th percentile** of the scaled counts.The value that 75% of scaled data are equal or less than that.
- Compute the average of these 75th percentiles across all samples.
- Select the reference sample whose 75th percentile is **closest to this average**.

This procedure helps identify a sample with a typical expression distribution.

---

### 3. Select genes for calculating the scaling factor  
*(Done separately for each sample relative to the reference)*

Different samples may use different sets of genes to derive their scaling factors. The goal is to use genes with the **least bias**, assuming that most genes are not differentially expressed.

#### a) Calculate log2 fold change (M value)
For each gene, calculate the **log2 fold change (M value)** between the sample and the reference.  
Genes with zero counts in either the sample or the reference are removed to avoid infinite values.

#### b) Calculate average expression (A value)
For each gene, calculate the **average expression (A value)** across the sample and the reference.  
This step helps identify genes with very high or very low expression.

#### c) Sort genes
Genes are sorted based on:
- **M values** (log fold changes)
- **A values** (average expression)

#### d) Apply trimming
- Remove genes with extreme M values (typically the top and bottom ~30%), which are likely biased or truly differentially expressed.
- Remove genes with extreme A values (typically the top and bottom ~5%), representing very highly or very lowly expressed genes.

#### e) Select remaining genes
The genes that remain after trimming are assumed to be unbiased and are used to estimate the scaling factor.

#### f) Repeat for each sample
This gene selection and trimming process is performed independently for each sample relative to the reference.

---

### 4. Calculate the weighted trimmed mean of log2 ratios

Using the remaining genes, edgeR calculates a **weighted trimmed mean of M values**, which is the core idea behind **TMM (Trimmed Mean of M values)**.

Genes with higher counts receive more weight because they provide more reliable information, while low-count genes contribute less.

---

### 5. Convert from log2 scale to normal scale

The weighted average of log2 ratios is converted back to the normal scale, producing the **normalization (scaling) factor** for that sample.

---

### 6. Repeat for all samples

This entire procedure is repeated for every sample relative to the reference sample.

---

### 7. Center scaling factors around 1

Finally, all scaling factors are divided by their **geometric mean**, ensuring that the normalization does not change the overall scale of the data across samples.

---

## Summary

edgeR’s TMM normalization corrects for both **library size differences** and **compositional bias** by estimating sample-specific scaling factors, which are later incorporated as offsets in the negative binomial model used for differential expression testing.

---

🔍 Watch this video by Josh Starmer, steps are explained clearly: 
https://www.youtube.com/watch?v=Wdt6jdi-NQo&list=PLblh5JKOoLUJo2Q6xK4tZElbIvAACEykp&index=12
