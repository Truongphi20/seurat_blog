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

The algorithm begins with the input feature matrix $A$ ($n$ cells $\times$ $m$ features) and an initial normalized vector $p_1 \in \mathbb{R}^m$ representing random feature residuals of an arbitrary cell.

```{math}
:label: lanczos-begin
\begin{cases}
\begin{aligned}

\alpha_1 = \lVert A \cdot p_1 \rVert \\
q_1 = \frac{1}{\alpha_1}A \cdot p_1  \\


\end{aligned}
\end{cases}
```

At the start of the iteration, $A \cdot p_1$ computes the cell similarity scores between all cells in $A$ and the arbitrary cell. Its magnitude is stored as $\alpha_1$, and the scaled unit vector $q_1 \in \mathbb{R}^n$ as described in Equation [](#lanczos-begin).

#### The loop runs

The number of iterations $k$ defines the size of the working Krylov subspace, as reflected by the number of desired singular vectors $v$ (by default, $v=50$; $k=v+7$).

```{math}
:label: lanczos-baseline
\bar{r}_j = A^T q_j - \alpha_j p_j
```

For iteration $j$, the un-orthogonalized feature residual vector $\bar{r}_j \in \mathbb{R}^m$ is defined as Equation [](#lanczos-baseline). Intuitively, $\bar{r}_j$ represents the updated feature-weighted covariate score ($A^T q_j$) minus the baseline covariate score ($\alpha_j p_j$).

```{math}
:label: lanczos-orthogon
r_j = \bar{r}_j - P_j (P_j^T \bar{r}_j)
```

Subsequently, Gram-Schmidt orthogonalization is performed on $\bar{r}_j$ as shown in Equation [](#lanczos-orthogon) (where $P_j = [p_1, p_2, \dots, p_j]$), extracting a perpendicular residual vector $r_j$ that is strictly independent from all previously computed feature residual vectors ($P_j$).

:::{tip} Why does $PP^Tv$ represent residual of $v$ on $P$-space? 

![](./static/perpendicular.png)

```{math}
\begin{aligned}

P^{T}e &= 0 \\
P^{T}(v - P\hat{x}) &= 0 \\
P^{T}P\hat{x} &= P^{T}v \\

\end{aligned}
```

```{math}
\begin{cases}
\begin{aligned}

\hat{x} &= P^{T}v \\
\hat{r} &= PP^{T}v

\end{aligned}
\end{cases}
```

:::

```{math}
\begin{cases}
\begin{aligned}

p_{j+1} &= \frac{r_j}{\lVert r_j \rVert} \\
\bar{q}_{j+1} &= A \cdot p_{j+1} - \lVert r_j  \rVert q_j

\end{aligned}
\end{cases}
```


```{math}
\hat{q}_{j+1} = \bar{q}_{j+1} - Q_{j+1} (Q^{T}_{j+1} \bar{q}_{j+1})
```

```{math}
\begin{cases}
\begin{aligned}

\alpha_{j+1} = \lVert \hat{q}_{j+1} \rVert \\
q_{j+1} = \frac{\hat{q}_{j+1}}{\alpha_{j+1}}

\end{aligned}
\end{cases}
```

#### Terminating

### Augmented Lanczos Bidiagonalization


## Summary