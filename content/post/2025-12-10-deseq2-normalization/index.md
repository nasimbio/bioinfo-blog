---
title: "How Does DESeq2 Perform Normalization?"
author: "Nasim Rahmatpour"
date: "2025-12-10"
slug: "deseq2-normalization-explained"
draft: false
categories: ["RNA-seq", "Normalization"]
tags: ["Normalization", "RNA-seq", "DESeq2", "DGE", "RNA-seq", "bulk-RNA"]
math: true
---

Differential gene expression (DGE) analysis requires correcting for two major
sources of technical variation:

1. **Library size differences** (sequencing depth)
2. **Compositional bias** between samples

---

## Why gene length normalization is *not* needed?

Gene length normalization (like normalization methods in RPKM or TPM) is **not required** for
DGE analysis because:

- DGE compares the **same gene across conditions**, not different genes within a sample and the same gene has the same length in different samples.

---

## What is compositional bias and why does it matter?

Compositional bias occurs when a subset of genes is **highly expressed in one
condition or tissue but not in others**.

For example:

- In one tissue, a small set of tissue-specific genes may dominate the read counts.
- Because sequencing depth is fixed, this forces other genes to receive fewer reads,
  even if their true expression has not changed.
- Simply scaling by total library size would make many genes appear falsely
  down-regulated.

DESeq2 normalization is specifically designed to correct for **both library size
differences and compositional bias**. To adjust for library size and compositional bias, it calculates a size factor.

---

## DESeq2 normalization: the median-of-ratios method

DESeq2 estimates a **size factor** for each sample using the
*median-of-ratios* approach. These size factors are later used to normalize raw counts.

---

## Step-by-step explanation

### 1. Log-transform raw counts

For each gene in each sample, DESeq2 considers the natural logarithm of the raw
read counts.  

---

### 2. Compute a pseudo-reference using the geometric mean

Then it calculates the average of logs for each gene across different samples. This is called **geometric mean**. This value represents a pseudo-reference expression level for that gene.

---

### 3. Exclude genes with zero counts

Genes that have Infinity after log transformation (or zero counts before log transformation) in one or more samples are excluded for downstram analysis.

---

### 4. Compute log ratios

For each gene *g* in each sample *s*, DESeq2 computes the log ratio:

`log(count_{g,s}) - log(geoMean_g)`

---

### 5. Take the median across genes

For each sample, DESeq2 takes the **median** of all gene-level log ratios.
This single value summarizes the overall shift of the sample relative
to the reference.

---

### 6. Convert to a size factor

Convert the median ratio to normal numbers. The median log ratio is exponentiated
to obtain the **size factor** for that sample:

`size_factor_s = exp(median_log_ratio_s)`

Here:

- `size_factor_s` is the size factor for sample *s*
- `median_log_ratio_s` is the median log ratio for sample *s*

---

### 7. Normalize the counts

Finally, raw counts are normalized by dividing by the sample-specific size factor:

`normalized_count_{g,s} = raw_count_{g,s} / size_factor_s`

These normalized counts are then used for downstream modeling and differential
expression testing.

---

## Why this approach works well

The median-of-ratios method is robust because:

- **Geometric means** reduce the influence of very highly expressed genes.
- **Medians** are resistant to outliers.
- Normalization is driven by **moderately expressed genes** like **housekeeping genes**, which better
  represent global expression behavior.

This makes DESeq2’s normalization particularly effective at correcting
compositional bias while preserving true biological differences.

---

🔍 For better visualization examples, see this link:  
*(https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)* 
