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

```{math}
:label: loess
\large \hat{\sigma_{i}}^2 = \text{LOESS} \left( \log_{10}(\sigma_i^2) \sim \log_{10}(\mu_i) \right)
```

Subsequently, the expected variance ($\hat{\sigma_{i}}^2$) is estimated by the Local Polynomial Regression model (LOESS) [@cleveland2017local]. The model fits a smooth trend to capture the relationship between gene abundance and variance in log-log space (the formula [](#loess)). Using local parabolic fitting, the expected variance is estimated across the mean expression range to generate a continuous, smooth curve (See [brilliant Josh's explanation](https://www.youtube.com/watch?v=Vf7oJ6z2LCc)).

```{math}
:label: compute-std-var
\Large
\begin{cases}
\begin{aligned}
z_{ij} &= \frac{c_{ij} - \mu_{i}}{\hat{\sigma_{i}}} \\
\bar{\sigma_{i}}^2 &= \frac{\sum_{j=1}^{N}{\left[\min(\sigma_{\text{max}}, z_{ij})\right]^2}}{N-1}
\end{aligned}
\end{cases}
```

Finally, the standardized variance ($\bar{\sigma_{i}}^2$) for gene $i$ is computed by formula [](#compute-std-var), which represents the variance of the standardized values ($z_{ij}$) across all cells, capped at a maximum value of $\sigma_{\text{max}} = \sqrt{N}$, which ensures the balanced distribution of variance, preventing outliers from overwhelmingly escalating the gene's overall variance.

The genes are then ranked by their standardized variance in descending order to select the top highly variable features (default `nfeatures = 2000`).

:::{tip}Why not utilize raw variance directly for sorting out features?

Due to the fact that variance rapidly increases as expression values rise [@ahlmann-eltzeComparisonTransformationsSinglecell2023], even when utilizing variance calculated from log-normalized data, highly expressed genes will dominate the top rankings [@stuartComprehensiveIntegrationSingleCell2019]. By standardizing the variance instead, this effect is controlled while retaining relative, feature-specific variability.

:::

### Mean variance plot (`mvp`)

```{math}
:label: mean-mvp
\Large \mu_{i} = ln \left( \frac{1}{N} \sum^{N}_{j=1}{e^{LN_{ij}}} \right)
```


```{math}
\Large
\begin{cases}
\begin{aligned}
  rm_i   &= \frac{1}{N} \sum_{j=1}^{N} \left( e^{LN_{ij}} - 1 \right) \\
  SS_i   &= \sum_{j=1}^{N} \left( e^{LN_{ij}} - 1 \right)^2 \\
  var_i  &= \frac{1}{N - 1} \left( SS_i - N \cdot rm_i^2 \right) \\
  disp_i &= \ln \left( \frac{var_i}{rm_i} \right)
\end{aligned}
\end{cases}
```


### Dispersion (`disp`)


## Summary