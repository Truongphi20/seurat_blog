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

Firstly, the variance ($\sigma_{i}^2$) corresponding to each gene is computed following the formula [](#compute-var), where $c_{ij}$ is an element of gene $i$ and cell $j$ in the count matrix; $\mu_{i}$ is the mean of expression count of gene $i$; and N is the number of cells.

Subsequently, the expected variance ($\hat{\sigma_{i}}^2$) is estimated by the Local Polynomial Regression model (LOESS) [@cleveland2017local]. The model fits a smooth trend to capture the relationship between gene abundance and variance in log-log space ($\log_{10}(\sigma_i^2) \sim \log_{10}(\mu_i)$). Once the polynomial model is fitted, the expected variance is calculated by passing the observed mean value into the model.

```{math}
:label: compute-std-var
{\Large
\begin{aligned}
z_{ij} &= \frac{c_{ij} - \mu_{i}}{\hat{\sigma_{i}}} \\
\bar{\sigma_{i}}^2 &= \frac{\sum_{j=1}^{N}{\left[\min(\sigma_{\text{max}}, z_{ij})\right]^2}}{N-1}
\end{aligned}
}
```

Finally, the standardized variance ($\bar{\sigma_{i}}^2$) for gene $i$ is computed by formula [](#compute-std-var), which represents the variance of the standardized values ($z_{ij}$) across all cells, capped at a maximum value of $\sigma_{\text{max}} = \sqrt{N}$. The genes are then ranked by their standardized variance in descending order to select the top highly variable features (default `nfeatures = 2000`).

### Mean variance plot (`mvp`)

### Dispersion (`disp`)


## Summary