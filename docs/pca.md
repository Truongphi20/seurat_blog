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
\bar{r}_{j+1} = A^T q_j - \alpha_j p_j
```

For iteration $j$, the un-orthogonalized feature residual vector $\bar{r}_{j+1} \in \mathbb{R}^m$ is defined as Equation [](#lanczos-baseline). Intuitively, $\bar{r}_j$ captures the updated feature loadings ($A^T q_j$) after stripping out the scaled baseline loadings ($\alpha_j p_j$) carried over from the previous feature loading.

```{math}
:label: lanczos-orthogon
r_{j+1} = \bar{r}_j - P_j (P_j^T \bar{r}_{j+1})
```

Subsequently, Gram-Schmidt orthogonalization is performed on $\bar{r}_{j+1}$ as shown in Equation [](#lanczos-orthogon) (where $P_j = [p_1, p_2, \dots, p_j]$) to ensure the next feature axis $p_{j+1}$ is strictly independent of all previous feature vectors in $P_j$.

:::{tip} Why does $PP^Tv$ represent residual of $v$ on the $P$-space? 

![](./static/perpendicular.png)

Consider an arbitrary vector $v$ and a linear vector dimension $P$, $v$ always includes 2 components: the residual $r$ on $P$, which is constructed from vectors in the $P$-space with coefficients $\hat{x}$; and the perpendicular component $e$.

```{math}
:label: penper-coeff
\begin{aligned}

P^{T}e &= 0 \\
P^{T}(v - P\hat{x}) &= 0 \\
P^{T}P\hat{x} &= P^{T}v \\

\end{aligned}
```

Hence $e$ is orthogonal (vertical) against $P$, the product of $e$ and foundational vectors of $P$ is 0. Based on that, the relationship between coefficents $\hat{x}$ and vector $v$ is described as Equation [](#penper-coeff), which is well known as the [Normal Equation](https://www.geeksforgeeks.org/machine-learning/ml-normal-equation-in-linear-regression/).    

```{math}
:label: penper-residual
\begin{cases}
\begin{aligned}

\hat{x} &= P^{T}v \\
r &= PP^{T}v

\end{aligned}
\end{cases}
```

Due to the presumption that $P$ is orthogonal, $P^{T}P = I$. Therefore, the coefficients $\hat{x}$ and residual $r$ on the $P$ dimension of vector $v$ are shown as Equation [](#penper-residual).

:::

```{math}
:label: q-raw-bar
\begin{cases}
\begin{aligned}

p_{j+1} &= \frac{r_{j+1}}{\lVert r_{j+1} \rVert} \\
\bar{q}_{j+1} &= A \cdot p_{j+1} - \lVert r_{j+1}  \rVert q_j

\end{aligned}
\end{cases}
```

Similar to Equation [](#lanczos-baseline), after the new feature axis $p_{j+1}$ is determined by normalizing feature-load residual $r_{j+1}$, the next unothogonalized cell residual $\bar{q}_{j+1}$ is computed from the cell loading ($A \cdot p_{j+1}$) being dependent from the previous estimated cell loading ($\lVert r_{j+1}  \rVert q_j$). The details is described as Equation [](#q-raw-bar).

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