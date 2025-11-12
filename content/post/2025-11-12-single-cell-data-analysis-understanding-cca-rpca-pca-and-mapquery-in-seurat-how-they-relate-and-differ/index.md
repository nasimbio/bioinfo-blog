---
title: "Single-Cell Data Analysis: Understanding CCA, RPCA, PCA, and MapQuery in Seurat; How They Relate and Differ"
author: "Nasim Rahmatpour"
date: "2025-11-12"
slug: "single-cell-data-analysis-understanding-cca-rpca-pca-and-mapquery-in-seurat-how-they-relate-and-differ"
draft: false
tags: ["single-cell", "Seurat", "CCA", "RPCA", "label transfer"]
categories: ["tutorial"]
---

**By Nasim**

🔹 1. CCA for Integration (Batch Effect Removal)
Canonical Correlation Analysis (CCA) identifies correlated gene expression patterns between datasets, finding the directions that maximize covariation across them.
Seurat uses these correlations to define anchors or pairs of biologically similar cells across datasets.
These anchors then drive the creation of a new integrated assay, which contains batch-corrected expression values suitable for joint analysis (clustering, visualization, and DGE).

🔹 2. CCA and PCA for Label Transfer
For label transfer, the goal shifts from merging datasets to mapping a new (query) dataset onto a labeled (reference) dataset.
Here, Seurat again finds anchors, but instead of correcting expression values, it uses them for KNN-based label propagation: each query cell inherits information (e.g., cell type, PCA coordinates) from its k nearest anchors in the reference space.
•	CCA-based label transfer: anchors are found via correlated features (CCA).
•	PCA-based label transfer: anchors are found using mutual nearest neighbors (MNNs) in PCA space, followed by the same KNN-weighted propagation.

🔹 3. RPCA for Faster Integration (Batch Effect Removal)
Reciprocal PCA (RPCA) is a more scalable, memory-efficient variant of CCA for integration.
Instead of computing a single joint CCA across all datasets, RPCA performs PCA separately for each dataset, then identifies MNN-based anchors across those PCA spaces.
This approach is ideal for large or closely related datasets.

🔹 4. MapQuery 
Unlike integration, MapQuery() is used after anchors are found to project query cells into the reference’s PCA or UMAP space and transfer labels.
It doesn’t merge datasets, instead, it aligns query data to the reference coordinate system.

Summary:
•	CCA → Finds correlated features between datasets; used for both integration and label transfer.

•	RPCA → PCA-based integration using MNN anchors; faster and scalable.

•	Integration → Produces an integrated assay with batch-corrected expression values.

•	Label transfer → Uses anchors + KNN-weighted propagation to map query data to a labeled reference.

•	CCA / RPCA → integration (new batch-corrected assay)

•	CCA / PCA + MapQuery → label transfer (projected query, label propagation)

