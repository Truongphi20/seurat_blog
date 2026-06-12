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

Based on the formula `~1`, a linear model matrix is assigned with only an intercept term for each sample. Additionally, a small ridge penalty ($\frac{10^{-10}}{N}$) is added to prevent overfitting.

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

Initially, a rough estimate of the dispersion for gene $i$ ($\widehat{\theta}_{i}$) is computed using Equation [](#rough-disper-est). Here, $\xi$ is an offset correction factor, defined as the reciprocal of the average exponentiated offset ($o_j$) across all $N$ samples. Since the offset in this context is the natural logarithm of the total counts for each sample, $\xi = (\overline{C})^{-1}$, where $\overline{C}$ is the mean total count across samples.

The rough dispersion estimator is derived from the mean–variance relationship of the Gamma-Poisson distribution ($\sigma^2 = \mu + \theta \mu^2$) [@ahlmann-eltzeGlmGamPoiFittingGammaPoisson2021], adjusted by the factor $\xi$ and constrained to be non-negative.

Here, $C_i$ denotes the vector of counts for gene $i$, while $var()$ and $mean()$ represent the variance and mean functions, respectively.

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
\beta_i^{(n+1)} = \beta_i^{(n)} + \frac{dl_{i}(\beta_i)}{dl'_{i}(\beta_i)}
```

Beta is estimated using the Newton-Raphson method [@akram2015newton], which iterates Equation [](#beta-est-step) until the value of beta converges. $dl_{i}$ is the derivative of likelihood function that count values of gene $i$ is of Negative Binomial (NB) distribution.  

```{math}
:label:beta-est-newton-raphson
\Large 
\begin{cases}
\begin{aligned}

\mu_{ij} &= \exp(\beta_i + o_j) \\
dl_{i}(\beta_i) &= \sum_{k=1}^{N}{\frac{C_{ik}-\mu_{ik}}{1 + \mu_{ik}  \theta_i}} \\
dl'_{i}(\beta_i) &= - \sum_{k=1}^{N}{\frac{\mu_{ik}(1+C_{ik}\theta_i)}{(1 + \mu_{ik}  \theta_i)^2}}

\end{aligned}
\end{cases}
```

Starting with an initial rough estimate of $\beta_i$ for gene $i$, the sample-specific expected count $\mu_{ij}$ is updated in each iteration $n$ by combining the biological parameter ($\beta_i$) and technical sampling offset ($o_j$). Even though $\mu$ represents the theoretical mean of the NB distribution, incorporating these sample-specific offsets allows the model to account for variation in sequencing depth across the samples. This helps strip away library-size biases while computing the log-likelihood function and its derivative according to Equation [](#beta-est-newton-raphson).

:::{tip} Maximum likelihood of NB distribution

By assuming counts of each gene follows NB distribution, where variance quadratically increases by mean, data is utillized to estimated coefficents in the model by the maximum likelihood method. 

```{math}
:label: nb-dis
\begin{cases}
\begin{aligned}

P(Y = y) &= \binom{y+r-1}{r-1} p^r (1-p)^y \\
p &= \frac{1}{1+\mu/r}

\end{aligned}
\end{cases}
```

Theoretically, the probability mass function (PMF) of NB distribution [](#nb-dis) estimates the probability of obtaining $y$ failures before archiving $r$ successes in a sequence of independent Bernoullian trials (two possible results, success or failure), with $p$ is the probability of succeeding any trial, and $\mu$ is the mean of number of failures [@sinharayDiscreteProbabilityDistributions2010]. In practical, since $y$ is a discrete value starts from 0, it is chosen to assign for $C_k$ without any biological meaning.     

```{math}
:label: likelihood-func
l(p,r) = \ln \left( \prod_{k=1}^{N}{ P(Y = C_k) } \right)

```

The log-likelihood function [](#likelihood-func) converts the product of probabilities into a sum, making it analytically straightforward to differentiate. As mentioned earlier, the relationship between the dispersion parameter $\theta$ and the number of successes $r$ is $r = 1/\theta$.

```{math} 
:label: likelihood-mu

l(\mu,\theta) = \sum_{k=1}^N{\left(  \ln \left[\binom{C_{k}+1/\theta-1}{1/\theta-1} \right] - \frac{1}{\theta}\ln(\mu\theta + 1) + C_k\ln(\frac{\mu\theta}{1+\mu\theta}) \right) }

```

Substituting $p$ according to [](#nb-dis), Equation [](#likelihood-mu) computes likelihood depending on the mean count of a gene. Maximum likelihood is found by rooting derivative of the likelihood function.  

```{math}
:label: likelihood-derivative
\frac{\partial l}{\partial \mu} = \sum_{k=1}^N{\frac{C_k - \mu}{\mu(1+\mu\theta)}}
```

Equation [](#likelihood-derivative) provides the partial derivative of the log-likelihood with respect to $\mu$.

```{math}
:label: likelihood-derivative-beta

\frac{\partial l}{\partial \beta} = \frac{\partial l}{\partial \mu} \cdot \frac{\partial \mu}{\partial \beta} = \sum_{k=1}^N{\frac{C_k - \mu}{1+\mu\theta}}

```

Hence the package `glmGamPoi` uses logarithm-scaled mean $\beta$ mentioned in [](#beta-est-newton-raphson) instead of count mean $\mu$, it is transformed to paritial derivative by $\beta$ following Equation [](#likelihood-derivative-beta).

:::

### Overdispersion estimation

```{math}
l(\theta) = \sum_{k=1}^{N}{ \left[ \ln\Gamma(C_k + \theta^{-1}) - \ln\Gamma(\theta^{-1}) - \left( C_k +\theta^{-1} \right)\ln(\mu_k+\theta^{-1}) - \frac{1}{\theta} \ln(\theta) \right] }
```

```{math}
\frac{d l}{d \theta} = \frac{1}{\theta} \left[ \sum_{k=1}^{N}{ \left( -\frac{1}{\theta} \left( \psi(C_k + \theta^{-1}) - \psi(\theta^{-1})  \right)  + \ln(1+\mu_k\theta) + \frac{C_k - \mu_k}{\mu_k + \theta^{-1}} \right) } \right]
```

```{math}
- \underbrace{\frac{1}{\theta}\sum_{k=1}^{N}{ \left( \psi(C_k + \theta^{-1}) - \psi(\theta^{-1}) \right) } }_{D(\theta)} + \underbrace{\sum_{k=1}^{N}{ \left( \ln(1+\mu_k\theta) + \frac{C_k - \mu_k}{\mu_k + \theta^{-1}} \right) } }_{L(\theta)} = 0
```

```{math}
G(\theta) = L(\theta) - D(\theta)
```

```{math}
\frac{d^{2}l}{d\theta^2} = -\theta^{-2}G(\theta) + \theta^{-1} G'(\theta)
```

:::{tip} The Cox-Reid (CR) adjustment
```{math}
\Large
\begin{cases}
\begin{aligned}

w_i &= \frac{1}{\frac{1}{\mu_i}+\theta} \\
B &= M^T \cdot IW \cdot M \\
B^{-1} &= (B + \lambda I)^{-1} \\
CR &= \frac{\partial}{\partial \theta} \left( -\frac{1}{2} \log |B| \right) = -\frac{1}{2}\mathrm{Tr}\left( B^{-1} \frac{\partial B}{\partial \theta} \right)

\end{aligned}
\end{cases}
```

```{math}
\Large
\begin{cases}
\begin{aligned}

DG = \sum_{y}{Fr(y) \cdot \psi(y + 1/\theta)}

\end{aligned}
\end{cases}
```
:::

### Shrinkage

## Summary