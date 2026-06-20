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

### Rough dispersion estimation

```{math}
:label:rough-disper-est
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
\beta_i^{(n+1)} = \beta_i^{(n)} + \frac{dl_{i}(\beta_i)}{dl'_{i}(\beta_i)}
```

Beta is estimated using the Newton-Raphson method [@akram2015newton], which iterates Equation [](#beta-est-step) until the value of beta converges. $dl_{i}$ is the derivative of likelihood function that count values of gene $i$ is of Negative Binomial (NB) distribution.  

```{math}
:label:beta-est-newton-raphson 
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

### Overdispersion (theta) estimation

Hence using the likelihood function $l(\mu,\theta)$ to model NB distribution for each gene, the overdispersion parameter $\theta$ is estimated after determination of $\beta$ derived from $\mu$. The vector of mean count $\mu_i$ of gene $i$ is deduced via the equation in [](#beta-est-newton-raphson), which stands for constant count mean across $N$ samples during optimizing $\theta$.      

```{math}
:label: likelihood-theta
l(\theta) = \sum_{k=1}^{N}{ \left[ \ln\Gamma(C_k + \theta^{-1}) - \ln\Gamma(\theta^{-1}) - \left( C_k +\theta^{-1} \right)\ln(\mu_k+\theta^{-1}) - \frac{1}{\theta} \ln(\theta) \right] }
```

The likelihood function $l(\theta)$ presented in Equation [](#likelihood-theta) is derived from general likelihood function of NB distribution [](#likelihood-mu). Here the factorial components in the combination term is substituted with Gamma functions to generalize the domain to continuous real numbers ($\mathbb{R}$), while any additive terms independent of $\theta$ are omitted.

```{math}
:label: likelihood-derivative-theta
\frac{d l}{d \theta} = \frac{1}{\theta} \left[ \sum_{k=1}^{N}{ \left( -\frac{1}{\theta} \left( \psi(C_k + \theta^{-1}) - \psi(\theta^{-1})  \right)  + \ln(1+\mu_k\theta) + \frac{C_k - \mu_k}{\mu_k + \theta^{-1}} \right) } \right]
```

To maximize the log-likelihood function, we first compute its derivative with respect to $\theta$, as shown in Equation [](#likelihood-derivative-theta), where $\psi(x)$ denotes the digamma function (the first derivative of $\ln\Gamma(x)$).   

```{math}
:label: likelihood-derivative-theta-simplify
- \underbrace{\frac{1}{\theta}\sum_{k=1}^{N}{ \left( \psi(C_k + \theta^{-1}) - \psi(\theta^{-1}) \right) } }_{D(\theta)} + \underbrace{\sum_{k=1}^{N}{ \left( \ln(1+\mu_k\theta) + \frac{C_k - \mu_k}{\mu_k + \theta^{-1}} \right) } }_{L(\theta)} = 0
```

Setting $dl/d\theta = 0$ yields the root-finding problem in Equation [](#likelihood-derivative-theta-simplify), which can be separated into two distinct components. The first part, $D(\theta)$, accounts for the digamma terms and involves intensive computation. To eliminate redundant calculations for identical count values, `glmGamPoi` optimizes this step by utilizing a count frequency table [@ahlmann-eltzeGlmGamPoiFittingGammaPoisson2021]. Meanwhile, the remaining log-mean component, $L(\theta)$, is more straightforward to compute.  

```{math}
:label: g-theta
G(\theta) = L(\theta) - D(\theta)
```

Structurally, the gradient $G(\theta)$ defined in Equation [](#g-theta) governs the direction of the variation of $\theta$. When $G(\theta) > 0$, the log-likelihood slope is positive, indicating that $\theta$ must be increased. Conversely, when $G(\theta) < 0$, the slope is negative, meaning $\theta$ must be decreased. A state of $G(\theta) = 0$ indicates that the log-likelihood has reached a local extremum.

To verify whether this stationary point corresponds to a local maximum or minimum. The second-derivative of the gradient is tested, if $d^{2}G/d\theta^2 < 0$, the curvature is concave down, confirming that the extremum is a local maximum representing the Maximum Likelihood Estimate (MLE) of $\theta$.

:::{caution}
The second derivative of log-likelihood is not compatible to the code implement, which is declared at https://github.com/const-ae/glmGamPoi/issues/74.
:::

:::{tip} The Cox-Reid (CR) adjustment

```{math}
\begin{cases}
\begin{aligned}

w_{jj} &= \frac{1}{\theta + 1/\mu_j } \\
B &= M^TWM \\
b(\theta) &= \ln(det(B))

\end{aligned}
\end{cases}
```

Based on the formula `~1`, a linear model matrix ($M$) is assigned with only an intercept term for each sample (a vector of 1). Additionally, a small ridge penalty ($\frac{10^{-10}}{N}$) is added to prevent overfitting.


```{math}

l(\theta)_{cr} = l(\theta) - \frac{1}{2} b(\theta) F_{cr}
```


:::

### Shrinkage

Even though NB distribution is able to obtain flexible variance-mean relationship, it is proved virtually being uncertainty in parameter estimation. Quasi-likelihood (QL) comes as a solution to adjust dispersion [@lundDetectingDifferentialExpression2012]. The variance-mean equation of quasi-distribution by genes is shown as [](#quasi-likelihood).    

```{math}
:label: quasi-likelihood
\sigma^2 = \theta_{QL}(\mu + \mu^2\theta_{\text{trend}})
```

Initially, the trend dispersion ($\theta_{\text{trend}}$) is derived from the vector of likelihood estimated dispersion ($\exp(1/\theta_{\text{ML}})$), which is performed local median regression to stablize dispersion. Specifically, the dispersion according to each mean count value ($\mu_i$) would be performed weighted median [@cormen2022introduction] with neighbors (sorted by mean, analyzing 100 neighbors by default), where the weight range originates from the range of probability in the domain $[-3,3]$ of standard normal distribution. 


```{math}
:label: theta-ql
\theta_{QL} = \frac{1+\mu\theta_{\text{ML}}}{1+\mu\theta_{\text{trend}}}
```

$\theta_{QL}$ is defined by [](#theta-ql), which derived from [](#quasi-likelihood) and mean-variance equation of NB distribution, where likelihood dispersion $\theta_{\text{ML}}$ was smoothed into $\theta_{\text{trend}}$. 

```{math}
:label: shrink_quasi_theta
\theta_{SQL} = \frac{\text{df}_0\tau^2_0 + \text{df}\theta_{QL}}{\text{df}_0 + \text{df}}
```

The shrunken quasi-likelihood overdispersion ($\theta_{SQL}$) for each gene is estimated by Equation [](#shrink_quasi_theta). Particularly, $\text{df}$ is the residual degree of freedom - DOF ($\text{\#samples} - 1$), while $\text{df}_0$ and $\tau^2_0$ are prior parameters representing DOF and dispersion scaler respectively. These prior parameters are estimated by modeling the $\theta_{QL}$ values using a scaled Chi-squared distribution via a natural cubic spline framework. The process successfully yeilds a theoretical pseudo-sample size ($\text{df}_0$) and a denoised central dispersion trend ($\tau_0^2$). Ultimately, the resulting $\theta_{SQL}$ reflects the true amount of informative dispersion captured from the data.  

## Summary