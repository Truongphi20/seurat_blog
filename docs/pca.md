---
title: Principal Component Analysis
numbering:
    math: true
---

## Introduction

:::{tip} Seurat command

```R
## Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

## Determine the ‘dimensionality’ of the dataset
ElbowPlot(pbmc)
```

:::

## Workflow

### Lanczos Bidiagonalization

#### Initializing

The algorithm begins with the input feature matrix $A$ ($n$ cells $\times$ $m$ features) and a random expression vector $p_1$ from an arbitrary cell.

```{math}
:label: lanczos-begin
\begin{cases}
\begin{aligned}

\alpha_1 = \lVert A \cdot p_1 \rVert \\
q_1 = \frac{1}{\alpha_1}A \cdot p_1  \\


\end{aligned}
\end{cases}
```

At the begin of iterations, $A \cdot p_1$ represents the covariance score of the arbitrary cell and cells in the $A$-base. Its magnitude $\alpha_1$ and scaled vector $q_1$ are described as Equation [](#lanczos-begin). 

#### The loop runs

```{math}
:label: lanczos-baseline
\bar{r}_j = A^T q_j - \alpha_j p_j
```

```{math}
:label: lanczos-orthogon
r_j = \bar{r}_j - p_{(1:j)} (p_{(1:j)}^T \bar{r}_j)
```

### Augmented Lanczos Bidiagonalization


## Summary