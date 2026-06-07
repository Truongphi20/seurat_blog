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
:label:beta-est-init
\Large 
\begin{cases}
\begin{aligned}

NC_{ij} &= \frac{C_{ij}}{\sum_{k=1}^{M}{C_{kj}}} \\
\widehat{\beta_i} &= \ln(mean(NC_i))

\end{aligned}
\end{cases}
```

At the beginning, rough beta of each gene ($\beta_i$) is estimated by Equation [](#beta-est-init), defined as the logarithm of the gene-oriented mean of normalised counts ($NC_{ij}$), which are proportion of count on total counts of each sample ($M$ is the number of genes).     

```{math}
:label:beta-est-step
\Large
\beta_i^{(n+1)} = \beta_i^{(n)} + \frac{dl_{i}}{ddl_{i}}
```

Beta is estimated using the Newton-Raphson method [@akram2015newton], which iterates Equation [](#beta-est-step) until the value of beta converges. $dl_{i}$ is the derivative of likelihood function that count values of gene $i$ is of Negative Binomial (NB) distribution.  

```{math}
:label:beta-est-newton-raphson
\Large 
\begin{cases}
\begin{aligned}

\mu_{ij} &= \exp(\beta_i + o_j) \\
dl_{i}(\beta_i) &= \sum_{k=1}^{N}{\frac{C_{ik}-\mu_{ik}}{1 + \mu_{ik}  \theta_i}} \\
ddl_{i}(\beta_i) &= \sum_{k=1}^{N}{\frac{\mu_{ik}(1+C_{ik}\theta_i)}{(1 + \mu_{ik}  \theta_i)^2}}

\end{aligned}
\end{cases}
```

:::{tip} Maximum likelihood of NB distribution 
:class: dropdown
:open: true

By assuming counts of each gene follows NB distribution, where variance parabolically increases by mean, data is utillized to estimated coefficents in the model by the maximum likelihood method. 

```{math}
:label: nb-dis
\begin{cases}
\begin{aligned}

P(Y = y) &= \binom{y+r-1}{r-1} p^r (1-p)^y \\
p &= \frac{1}{1+\mu/r)}

\end{aligned}
\end{cases}
```

Theoretically, the probability mass function (PMF) of NB distribution [](#nb-dis) estimates the probability of obtaining $y$ failures before archiving $r$ successes in a sequence of independent Bernoullian trials (two possible results, success or failure), with $p$ is the probability of succeeding any trial, and $\mu$ is the mean of number of failures [@sinharayDiscreteProbabilityDistributions2010]. In practical, since $y$ is a discrete value starts from 0, it is chosen to assign for $C_k$ without any biological meaning.     

```{math}
:label: likelihood-func
l(p) = \ln \left( \prod_{k=1}^{N}{ P(Y = C_k) } \right)

```

Likelihood function [](#likelihood-func) is logarithm of PMF to convert products into manageable sums, making it analytically straightforward to differentiate. Particularly, count number ($C_k$) represents the number of failures of biological captures. As the mean–variance relationship mentioned ealier, $r$ is defined by $\theta$ ($r = 1/\theta$) when expanning the likelihood function. 

```{math} 
:label: likelihood-mu

l(\mu) = \sum_{k=1}^N{\left(  \ln \left[\binom{C_{k}+1/\theta-1}{1/\theta-1} \right] - \frac{1}{\theta}\ln(\mu\theta + 1) + C_k\ln(\frac{\mu\theta}{1+\mu\theta}) \right) }

```

Substituting $p$ according to [](#nb-dis), Equation [](#likelihood-mu) computes likelihood depending on the mean count of a gene. Maximum likelihood is found by rooting derivative of the likelihood function.  

```{math}
:label: likelihood-derivative
\frac{\partial l}{\partial \mu} = \sum_{k=1}^N{\frac{C_k - \mu}{\mu(1+\mu\theta)}}
```

Due to likelihood function here is used to estimated $\mu$, the partial derivative of likelihood by $\mu$ is performed as [](#likelihood-derivative). 

```{math}
:label: likelihood-derivative-beta

\frac{\partial l}{\partial \beta} = \frac{\partial l}{\partial \mu} \cdot \frac{\partial \mu}{\partial \beta} = \sum_{k=1}^N{\frac{C_k - \mu}{(1+\mu\theta)}}

```

Hence the package `glmGamPoi` uses logarithm-scaled mean $\beta$ mentioned in [](#beta-est-newton-raphson) instead of logarithmic mean $\mu$, it is transformed to paritial derivative by $\beta$ following Equation [](#likelihood-derivative-beta).

:::

### Overdispersion refining 

## Summary