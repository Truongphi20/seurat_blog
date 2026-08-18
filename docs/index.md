---
title: "Seurat Under the Hood"
description: "Deep dives into Seurat internals"
---

## Motivation

[Seurat](https://satijalab.org/seurat/) is a standard bioinformatic tool used to explore and analyze single-cell RNA sequencing (scRNA-seq) data. Its ["Guided Clustering Tutorial"](https://satijalab.org/seurat/articles/pbmc3k_tutorial) serves as a foundational resource for newcomers learning standard single-cell workflows.

While the tutorial provides a solid scaffold by demonstrating key commands and high-level concepts, it lacks deep explanations of the underlying statistical and algorithmic methods. This creates a gap for readers who want to understand what happens "under the hood".

This blog serves as an extension of the Seurat vignette, detailing the processes running behind these straightforward commands and explaining the reasoning behind each step.


## Content

The content in this blog strictly follows typical functional steps in the Seurat (v5.5.0) vignette, breaking down the main command and underlying mechanisms of each step. The method for knitting content of this blog detailed in [Methodology](./method.md). 

```{image} ./static/Overview.png
:alt: Overview
:width: 80%
:align: center
```

First, the raw count matrix (such as from 10X Genomics) undergoes pre-processing to reduce technical noise and unwanted variance. This traditional pipeline includes three key steps: [Normalization](./normalization.md), [Feature selection](./variable_features.md), and [Scaling data](./scaling.md). Alternatively, [SCTransform](./sctransform.md) can be performed as an all-in-one approach.

Next, purified matrix is performed Principal Component Analysis (PCA) to capture major axes of biological variation across cells, which facilitates cell clustering. From all of that, Uniform Manifold Approximation and Projection (UMAP) is employed to fine-tuning clusters.  

Finally, marker genes are identified each cells proving bases for cell-type determination.

::: {attention}

Sections on dimensionality reduction and marker gene labeling will be available soon!

:::