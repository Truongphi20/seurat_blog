---
title: Gamma-Poisson GLM
numbering:
    math: true
---

## Introduction

:::{tip} Seurat command
:class: dropdown
:open: true
This command is captured from the process of [SCTransform](./sctransform.md).
```R
fit <- glmGamPoi::glm_gp(data = umi,
                           design = as.formula(new_formula),
                           col_data = data,
                           offset = log_umi,
                           size_factors = FALSE)

fit$theta <- 1 / fit$overdispersions
```

Command explanation:

- `data`: The count matrix, genes $\times$ samples.
- `design`: The statistical model that fitted data. The model `~1`, which treats all samples are in a same group, is used in this example.
- `col_data`: The metadata of samples.
- `offset`: Additional constants of each genes to adjust size factors. In this example, it is natural log of the total counts of each sample.
- `size_factors`:  Scaling factors to normalise counts across samples. In the commands, size factors are not applied. 

:::

## Workflow

Base on the formula `~1`, a linear-model matrix is assigned with only an intercept for each sample. Additionally, small ridge penalties ($\frac{10^{-10}}{N}$) is added to avoid overfitting.

### Rough dispersion estimation

```{math}
:label:rough-disper-est
\Large 
\begin{cases}
\begin{aligned}
\xi &= \left( \frac{1}{N}\sum_{j=1}^{N}{\exp(o_j)} \right)^{-1} \\
\widehat{\theta}_{i} &= \max \left(\frac{var(C_i) - \xi mean(C_i)}{mean(C_i)^2}, 0 \right)
\end{aligned}
\end{cases}
```

Firstly, a rough estimate of the dispersion level for gene $i$ ($\widehat{\theta}_{i}$) is computed according to Equation [](#rough-disper-est), where $\xi$ is an offset correction factor, which defined as the reciprocal of the average exponential of the offset ($o_j$) on total number of samples ($N$). Hence the offset in this command is natural log of the total counts of each sample, so $\xi = (\overline{C})^{-1}$, where $\overline{C}$ is the mean of total count across samples.

The rough dispersion estimator is derived from the mean–variance relationship of the Gamma-Poisson distribution ($\sigma^2 = \mu + \theta \mu^2$) [@ahlmann-eltzeGlmGamPoiFittingGammaPoisson2021], adjusted by the offset factor $\xi$ and ensured to be positive.

Here, $C_i$ denotes the vector of count values for gene $i$, while $var()$ and $mean()$ represent the variance and mean functions, respectively.

:::{tip} Why is $\xi$ necessary?
This question hasn't known yet.
:::

### Beta estimation

```{math}
:label:beta-est
\Large 
\begin{cases}
\begin{aligned}

NC_{ij} &= \frac{C_{ij}}{\sum_{k=1}^{M}{C_{kj}}} \\
\widehat{\beta_i} &= \ln(mean(NC_i))

\end{aligned}
\end{cases}
```

### Overdispersion refining 

## Summary