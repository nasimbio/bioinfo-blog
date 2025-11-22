---
title: "Centering and Scaling Before PCA in Single Cell RNA-seq; Why It Matters"
author: "Nasim Rahmatpour"
date: "2025-11-13"
slug: "centering-and-scaling-before-pca-why-it-matters"
draft: false
tags: ["single-cell", "Seurat", "Centering", "Scaling", "PCA", "heatmap"]
categories: ["tutorial"]
---


🧠 Centering and Scaling Before PCA; Why It Matters

When performing PCA on gene expression data, centering and scaling aren’t just optional, they’re essential to ensure PCA captures true biological variation rather than numerical artifacts. PCA is powerful, but only when the data are centered and scaled correctly.
Otherwise, it captures the wrong kind of “variance.” 

🎯 Centering

PCA ideally looks for directions of maximum variance around the mean of each gene.
But if we don’t center the data, PCA instead measures variance around the origin (0), not around each gene’s mean.

This subtle difference matters a lot:

That means housekeeping genes which are highly expressed everywhere even if they barely vary appear to contribute huge variance simply because their mean values are far from zero dominate the first PC, even though they don’t separate cell types.
The result: the first principal component becomes dominated by highly expressed, low-variance genes, masking real biological structure.

🧩 High mean ≠ informative variance

✅ Centering removes mean bias so PCA focuses on true biological variation.

⚖️ Scaling

If we skip scaling, genes with large numeric variance (often noisy or highly expressed) overpower the analysis. Z-scoring gives every gene equal influence, so PCA highlights correlated expression patterns not raw magnitude.

🧩 Scaling removes numeric bias, not biological signal.

🔬 In Seurat

1️⃣ FindVariableFeatures() → picks biologically variable genes

2️⃣ ScaleData() → z-score (center + scale)

3️⃣ RunPCA() → captures real structure across cells

🎨 Same reason we z-score heatmaps:

To compare genes relative to their own mean and variance.

DoHeatmap() in Seurat already uses scaled data.

Summary:

PCA without centering/scaling shows who’s loudest.
PCA with centering/scaling shows who’s different.


🧬 Bulk RNA-seq vs. Single-cell RNA-seq: Why Scaling Decisions Differ

Whether you scale before PCA depends on type of data. In bulk RNA each sample already represent an average over many cells and variance across samples reflect biological conditions. So, for bulk RNA you usually don’t scale before PCA since absolute variance reflects biological meaning and it corresponds to PCA based on covariance. For single cell, you should scale before PCA otherwise PCA will be driven by highly expressed genes rather than meaningful biological expression. This corresponds to PCA based on correlation matrix. 

To learn more, check out this excellent explanation on “chatomics”, YouTube channel. I always learn from him.

https://www.youtube.com/watch?v=P7kj0GLTgS4
