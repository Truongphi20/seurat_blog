# Scaling data

## Introduction

As the tutorial mentioned, scaling data is a prerequisite step for dimensionality reduction like Principal Component Analysis (PCA), where data per gene across cells is centered and standardized. 

This method ensures that all genes are given equal weight and are placed on the same scale (z-scores), allowing for accurate downstream comparisons. 

:::::{tip} Seurat command

```R
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

:::::

## Methods

The scaled-row value ($SR_{ij}$) is calculated for each gene $i$ and cell $j$. The data is scaled vertically by gene across all cells, using values from the normalized matrix ([log-normalize matrix](#log-normalize) in the case of the tutorial).

```{math}
:label: scaling-func
\Large
\begin{cases}
\begin{aligned}

\mu_i &= \frac{1}{N}\sum_{k=1}^{N}{LN_{ik}} \\
\sigma_i &= \sqrt{\frac{\sum_{k=1}^{N}{ (LN_{ik} - \mu_i)^2 }}{N-1}} \\
SR_{ij} &= min \left(\frac{LN_{ij} - \mu_i}{\sigma_i}, SR_{\text{max}} \right) 

\end{aligned}
\end{cases}
```

First, the mean value of each gene ($\mu_i$) is computed as the average of log-normalized values across all cells for gene $i$, where $N$ is the total number of cells. From that, standard deviation $\sigma_i$ is calculated for each gene. The normalized values are standardized using the mean and the standard deviation of each gene and clipped at a maximum threshold ($SR_{\max}$), which defaults to 10, controlled by the `scale.max` argument.        

## Summary

While [](./normalization.md) minimizes technical noise to enable accurate cell-to-cell comparisons, scaling standardizes the data to allow for direct comparisons across different genes. See [ScaleData.R](https://raw.githubusercontent.com/Truongphi20/seurat_blog/refs/heads/main/algorithm_code/ScaleData.R) for more details.