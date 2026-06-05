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
- `offset`: Additional constants of each genes to adjust size factors. In this example, it is natural log of the total counts of each gene.
- `size_factors`:  Scaling factors to normalise counts across samples. In the commands, size factors are not applied. 

:::

## Workflow

Base on the formula `~ 1`, a linear-model matrix is assigned with only an intercept for each sample. Additionally, small ridge penalties ($\frac{10^{-10}}{N}$) is added to avoid overfitting.

### Rough dispersion estimation

```{math}
:label:rough-disper-est
\Large RDE_{i} = \frac{var(C_i) - mean(C_i)/\overline{C}}{mean(C_i)^2}
```

Firstly, a rough estimate of the dispersion level for gene $i$ ($RDE_{i}$) is computed according to Equation [](#rough-disper-est), where $C_i$ represents the vector of raw count values for gene $i$, and $\overline{C}$ is the mean of the total counts across all samples. Note that $var()$ and $mean()$ denote the variance and mean functions, respectively.

### Beta estimation

### Overdispersion refining 

## Summary