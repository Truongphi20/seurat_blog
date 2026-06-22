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

At the begining, a number of genes is randomly selected (default is 2000) depending on quantity of expression. They are modeled by Gamma-Poisson general linear model (Gamma-Poisson GLM) to properly obtain overdispersion. See [](./glmGamPoi.md) for more detail about the fitting model process. 

Briefly, it uses the Maximum Likelihood Estimates (MLE) method to estimate coefficients of the model according to the count values of each gene. The coefficients includes relative log-scaled count mean $\beta$, and overdispersion $\theta$.

```{math}
\sigma^2 = \mu + \theta\mu^2
```

### Regularizing model

### Pearson residuals

## Summary