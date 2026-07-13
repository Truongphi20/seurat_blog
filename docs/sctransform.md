---
title: SCTransform workflow
numbering:
    math: true
---

## Introduction

:::{tip} Seurat command

```R
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
```

:::

## Workflow

### Fitting model

At the beginning, only genes obtains overdispersion factor ($\sigma^2 > \mu$) are used to train. A number of genes is randomly selected (default is 2,000) across their expression levels, and the total count of selected genes must be greater than 5 by default. Subsequently, their raw count matrix are modeled by Gamma-Poisson general linear model (Gamma-Poisson GLM) to properly obtain overdispersion. See [](./glmGamPoi.md) for more detail about the model-fitting process. 

Briefly, it uses the Maximum Likelihood Estimates (MLE) method to estimate coefficients of the model according to the count values of each gene. The coefficients include relative log-scaled count mean $\beta$, and overdispersion $\theta$. Note that the definition of overdispersion between sctransform and Gamma-Poisson GLM are reciprocal.

```{math}
:label: mean-var
\begin{cases}
\begin{aligned}

\sigma^2 &= \mu + \frac{1}{\hat{\theta}}\mu^2 \\
\alpha &= \hat{\theta} / \theta


\end{aligned}
\end{cases}
```

Next, theoretical overdispersion $\hat{\theta}$, which is derived from mean-variance relationship of NB model, shown as Equation [](#mean-var). $\alpha$ is the ratio of expected overdispersion ($\hat{\theta}$) over observed overdispersion ($\theta$), which originates from estimation above. If $\alpha < 0.001$, the model for that gene is assumed to follow the Poisson distribution ($\sigma^2 = \mu$). 

### Regularizing model

#### Geometric mean

```{math}
:label: geometric-mean-fn

\mu_{\text{g}} = \exp\left[ \frac{1}{N} \sum_{j=1}^{N}{\ln(C_j + \epsilon)} \right] - \epsilon

```

For each gene, logarithmic geometric mean ($\mu_{\text{g}}$) computed by Equation [](#geometric-mean-fn), represents as a centric count without affected by outlier samples. With $C_j$ is the count value of sample $j$ in total $N$ samples, and $\epsilon$ (default is 1) is a small fixed number to avoid $\ln(0)$ [@hafemeisterNormalizationVarianceStabilization2019].

#### Outlier detection

The local outlier genes followed the $\mu_{g}$-axis are detected based on the matrices of the log-scaled overdispersion factor ($F = 1+\mu/\theta$) and $\beta$ estimated from Gamma-Poisson GLM. They represent for outlier detectors of variance and mean perspectively.

Initially, the genes are isolated into bins according to their $\mu_{g}$ value with the heuristic binwidth rescaled from the optimal bandwidth by the method of @sheatherReliableDataBasedBandwidth1991 to ensure data fairly staying grounded on bins (https://github.com/satijalab/sctransform/issues/214).   

```{math}
:label: outlier-score

S = \frac{X - \text{median}(X)}{\text{MAD}(X)}

```

For each evaluated matrix ($X$ is $F$ or $\beta$), the outlier scores ($S$) are computed for the points in each bin following Equation [](#outlier-score), which measures the distance to the median and eliminates deviation by [Median Absolute Deviation (MAD)](https://en.wikipedia.org/wiki/Median_absolute_deviation) to enable comparability. For enhancing reliability, the measurement is processed on two grids lating on each other a half of binwidth. A outlier gene is labeled on a matrix when the scores on both grids are larger than a threshold, which defaults 10.

Ultimately, a gene is marked as an outlier gene when both evaluated matrices (variance and mean) agree.

### Pearson residuals

## Summary