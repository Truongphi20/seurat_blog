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

```{math}
:label: lanczos-begin
\begin{cases}
\begin{aligned}

\alpha_1 = \lVert Ap_1 \rVert \\
q_1 = \frac{1}{\alpha_1}Ap_1  \\


\end{aligned}
\end{cases}
```

```{math}
:label: lanczos-baseline
r_j = A^T q_j - \alpha_j p_j
```

```{math}
:label: lanczos-orthogon
r_j = r_j - p_{(1:j)} (p_{(1:j)}^T r_j)
```

### Augmented Lanczos Bidiagonalization


## Summary