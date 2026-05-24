---
numbering:
    math: true
---

# Feature selection

## Introduction

## Methods

### Variance Stabilizing Transformation (`vst`)

This method computes standardized variance based on the layer input (the `count` matrix according the tutorial).

```{math}
:label: compute-var
\Large \sigma_{i}^2 = \frac{\sum_{j=1}^{N}{(c_{ij} - \mu_{i})^2}}{N-1}
```

Firsly, the variance ($\sigma_{i}^2$) corresponding for each gene is computed following the formula [](#compute-var) with $c_{ij}$ is a element of gene $i$ and cell $j$ in the count matrix; $\mu_{i}$ is the mean of expression count of gene $i$; and N is the number of cell.

Subsequently, the expected variance ($\hat{\sigma_{i}}^2$) is estimated by the Local Polynomial Regression model [@cleveland2017local], which applies linear regression to determine likelihood variance on polinomial graph in a window along the mean value axis. When polynomial model is identified, expected variance is returned by the mean values.

```{math}
:label: compute-std-var
{\Large
\begin{aligned}
z_{ij} &= \frac{c_{ij} - \mu_{i}}{\hat{\sigma_{i}}} \\
\bar{\sigma_{i}}^2 &= \frac{\sum_{j=1}^{N}{\left[\min(\sigma_{\text{max}}, z_{ij})\right]^2}}{N-1}
\end{aligned}
}
```

Finally, the standardized variance ($\bar{\sigma_{i}}^2$) for gene $i$ is computed by formula [](#compute-std-var), which is the sum of square of standardized values ($z_{ij}$) across cells clipped by $\sigma_{\text{max}} = \sqrt{N}$. The genes are ranked by standardized variance decreasingly and getting the top 2000 most variable genes (default `nselect = 2000`). 

### Mean variance plot (`mvp`)

### Dispersion (`disp`)


## Summary