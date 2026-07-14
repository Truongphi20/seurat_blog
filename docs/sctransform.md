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

#### Local outlier detection

Local outlier genes along the $\mu_{g}$-axis are detected based on the matrices of the log-scaled overdispersion factor ($F = 1+\mu/\theta$) and $\beta$ estimated from Gamma-Poisson GLM. These parameters represent outlier indicators for variance and mean, respectively.

Initially, the genes are assigned to bins according to their $\mu_{g}$ value. The binwidth is determined using a heuristic rescale from the optimal bandwidth featured by @sheatherReliableDataBasedBandwidth1991 to ensure data is fairly distributed across bins (https://github.com/satijalab/sctransform/issues/214). 

```{math}
:label: outlier-score

S = \frac{X - \text{median}(X)}{\text{MAD}(X) + \epsilon}

```

For each evaluated matrix ($X$ representing either $F$ or $\beta$), the outlier scores ($S$) are computed for the points within each bin following Equation [](#outlier-score), which measures the distance to the median and eliminates deviation by [Median Absolute Deviation (MAD)](https://en.wikipedia.org/wiki/Median_absolute_deviation) to enable comparability. For enhancing reliability, the measurement is processed on two overlapping grids, offset by half a binwidth. A gene is flagged as an outlier within a specific matrix when its absolute score on both grids exceeds a threshold, which defaults to 10.

Ultimately, a gene is marked as an outlier if it is flagged by either of the evaluated matrices (variance or mean).

#### Kernel smoothing

```{math}
:label: smoothing-func
\begin{cases}
\begin{aligned}

\large w_{ij} &= K\left(\frac{|x_i - x_j|}{h}\right) \\
\large \overline{y}_j &= \frac{\sum_{i=i_{\text{min}}}^N y_i \cdot w_{ij} }{\sum_{i=i_{\text{min}}}^N w_{ij}} 

\end{aligned}
\end{cases}
```

### Pearson residuals

## Summary