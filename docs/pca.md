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

The algorithm begins with the input feature matrix $A$ ($n$ cells $\times$ $m$ features) and a random expression vector $p_1$ of an arbitrary cell.

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

The number of iteration for biagonalization $k$ is the size of working Krylov subspace reflected by the number of singular vectors to estimate $v$ ($v=50$ by default, $k=v+7$).

```{math}
:label: lanczos-baseline
\bar{r}_j = A^T q_j - \alpha_j p_j
```

For the iterator $j$, the based residual vector $\bar{r}_j$ defined as Equation [](#lanczos-baseline). Intuitively, it is the subtraction of the feature-specific covariate score and baseline covariate score of the based cell with the arbitrary cell.   

```{math}
:label: lanczos-orthogon
r_j = \bar{r}_j - P_{(1:j)} (P_{(1:j)}^T \bar{r}_j)
```

### Augmented Lanczos Bidiagonalization


## Summary