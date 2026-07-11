---
numbering:
    math: true
---

# SJ bandwidth estimation

## Introduction

[Kernel Density Estimation (KDE)](https://en.wikipedia.org/wiki/Kernel_density_estimation) is a common estimator allowing to model random-variance data, which does not follow any model assumptions (e.g. Bell curve of Normal Distribution). 

KDE mainly depends on the kernel function (usually Normal Distribution), which is assumed as local model for each point to calculate partitional probability for others based on their distances. Moreover, a smoothing parameter $h$ - bandwidth, is a important constant affecting to sensitivity of the model. If $h$ is too small, the model is oversmoothing (underfitting). In contrast, the estimator would noisily wiggle (overfitting) when the $h$ being too large [@Scott2012].

Eventhough many methods were developped to optimize $h$ value, the method of @sheatherReliableDataBasedBandwidth1991 is still a reliable approach by minimizing the erroneous differences of estimator for various complex nonparametric models [@eidousComparativeStudyBandwidth2010].     

:::{tip} Seurat command
This command is captured from the process of [SCTransform](./sctransform.md).

```R
bw.SJ(x)
```

Here `x` is an arbitrary numeric vector.
:::

## Workflow

Initially, the input data vector $X$ is partitioned into discrete bins to calculate a bin-step density (detailed in [](#bin-step-dens)). This binning strategy dramatically reduces the computational cost of evaluating the double summations required by the pilot bandwidth estimators when performing the secondary bandwidth $\alpha$ estimator detailed in [the secondary bandwidth ($\alpha$) estimator](#sec-bw-est). 

From that, the roughness functionals are estimated, allowing the optimal primary bandwidth $h$ to be found (as described in [the primary bandwidth $h$](#prim-bw-h-comp)) by solving the derivative of the AMISE equation, the core metric evaluating the accuracy of estimation. The underlying mathematical theories are detailed in the boxes below.

:::{tip} Why don't we just solve the minimization of AMISE to find bandwidth $h$?

Asymptotic mean integrated squared error (AMISE) measures dissimilarity between the underlying true density function and the estimator constructed from observed data. 

```{math}
:label: amise-eq
\text{AMISE}(h) = \underbrace{(nh)^{-1} R(K) \vphantom{\frac{1}{4}}}_{\text{variance}} + \underbrace{\frac{1}{4} h^4\sigma^4_K R(f'')}_{\text{bias}}
```

As detailed by @sheatherReliableDataBasedBandwidth1991, AMISE consists of two components representing the variance of the estimation and its squared bias relative to the true model:

- $n$: The sample size
- $h$: The smoothing bandwidth
- $R(g)$: The roughness functional ($\int{g(x)^2dx}$)
- $K$: The primary kernel function (typically a standard normal density)
- $f$: The true underlying density function

The ultimate target of the SJ method is choosing $h$ to minimize AMISE. While the variance term is mathematically straightforward, the true roughness $R(f'')$ in the bias term obstructs direct minimization because it depends on the unknown true density $f$. Crucially, $f$ is completely independent of our chosen bandwidth $h$.

```{math}
:label: amise-dev
h_{\text{opt}} = \left[ \frac{R(K)}{\sigma_K^4 R(f'') n} \right]^{1/5}
```

Therefore, the core workflow prioritizes to estimate $R(f'')$ first, allowing to find $h_{\text{opt}}$ by solving the minimized AMISE equation (Equation [](#amise-dev)) derived from setting its derivative to zero.

:::

:::{tip} How to estimate $R(f'')$?

The kernel-based method of @jonesUsingNonstochasticTerms1991 is utilized to estimate $R(f'')$.

```{math}
:label: roughness-define
R(f'') = \int_{-\infty}^{+\infty} [f''(x)]^2 dx = \int_{-\infty}^{+\infty} f^{\text{iv}}(x)f(x) dx = \mathbb{E}[f^{\text{iv}}(X)]
```

First, the roughness function of $f''$ is rewritten using integration by parts to express it as the mathematical expectation of the 4-th derivative of the true density $f$, as shown in Equation [](#roughness-define).

```{math}
:label: real-f-est
\hat{f}_\alpha^{\text{iv}}(x) = \frac{1}{n\alpha^5} \sum_{j=1}^n L^{\text{iv}}\left(\frac{x - X_j}{\alpha}\right)
```

To estimate this expectation, a secondary pilot kernel density estimator is introduced using a pilot kernel function $L$ and a pilot bandwidth $\alpha$ (distinct from the primary kernel $K$ and bandwidth $h$). Evaluating the 4-th derivative of this pilot density yields Equation [](#real-f-est).

```{math}
:label: s-diagonal
\hat{S}_{\text{D}}(\alpha) = \frac{1}{n} \sum_{i=1}^n \hat{f}_\alpha^{\text{iv}}(X_i) = \frac{1}{n^2 \alpha^5} \sum_{i=1}^n \sum_{j=1}^n L^{\text{iv}}\left(\frac{X_i - X_j}{\alpha}\right)
```

Evaluating this derivative estimator at every observed data point $X_i$ and taking the sample average gives the "diagonals-in" functional estimator $\hat{S}_{\text{D}}(\alpha)$ shown in Equation [](#s-diagonal).

```{math}
:label: sep-s-diag
\hat{S}_{\text{D}}(\alpha) = \frac{L^{\text{iv}}(0)}{n \alpha^5} + \hat{S}_{\text{ND}}(\alpha)
```

Separating $\hat{S}_{\text{D}}(\alpha)$ into two parts isolates the diagonal terms (where $i=j$) from the off-diagonal terms $\hat{S}_{\text{ND}}(\alpha)$ ($i \neq j$), which facilitates an asymptotic expansion of its expected value.

```{math}
:label: exp-s-diag
\mathbb{E}[\hat{S}_{\text{D}}(\alpha)] \approx R(f'') + \underbrace{\frac{L^{\text{iv}}(0)}{n \alpha^5}}_{\text{Positive Bias}} - \underbrace{\frac{1}{2}\alpha^2 \sigma_L^2 R(f''')}_{\text{Negative Bias}}
```

Equation [](#exp-s-diag) takes the expectation of both terms to obtain $\mathbb{E}[\hat{S}_{\text{D}}(\alpha)]$, which is the ultimate estimation for $\mathbb{E}[f^{\text{iv}}(X)]$. After applying the asymptotic expansion of $\mathbb{E}[\hat{S}_{\text{ND}}(\alpha)]$, there are approximally two components causing bias from the true roughness $R(f'')$. 


```{math}
:label: alpha-eq
\alpha = \left( \frac{2 L^{\text{iv}}(0)}{\sigma_L^2} \right)^{1/7} R^{-1/7}(f''') n^{-1/7}
```

To correct the estimator, the bias components is canceled out by the bandwidth $\alpha$ following Equation [](#alpha-eq), where positive and negative bias are equal each other. 

```{math}
:label: alpha-trans1
\alpha(h_{\text{opt}}) = \left( \frac{2 L^{\text{iv}}(0) \sigma_K^4}{\sigma_L^2 R(K)} \right)^{1/7} \cdot \left( \frac{R(f'')}{R(f''')} \right)^{1/7} h_{\text{opt}}^{5/7}
```

By manipulating with Equation [](#amise-dev), $\alpha$ is able to performed as Equation [](#alpha-trans1), where the ratio $R(f'')/R(f''')$ can be estimated by heuristic pilot bandwidths.

Summarily, $R(f'') \approx \hat{S}_{\text{D}}(\alpha(h_{\text{opt}}))$, which is performed as Equation [](#s-diagonal) with an independent KDE obtaining $\alpha$ depended on the optimized bandwidth $h_{\text{opt}}$ shown as Equation [](#alpha-trans1). From that, $h_{\text{opt}}$ can be solved by the equation of subsituting $\hat{S}_{\text{D}}(\alpha(h_{\text{opt}}))$ to Equation [](#amise-dev). 

:::

(bin-step-dens)=
### Bin-step density

First, the data is partitioned into $N$ discrete bins (the default is 1000) to determine the individual bin frequencies ($b_k$, where $0 \le k \le N-1$). These frequencies are used to efficiently compute the total number of data point pairs at each specific bin-step.

```{math}
:label: ci-sum

c_i = \sum_{k=i}^{N-1}{b_{k}b_{(k-i)}}
``` 

Equation [](#ci-sum) measures total density of unique pairs ($c_i$) between 2 bins with the $i$-bin distance. 

```{math}
:label: c0-sum
c_0 = \sum_{k=0}^{N-1}{ b_k \choose 2 } = \sum_{k=0}^{N-1}{\frac{b_k (b_k-1)}{2}}
```

For the baseline case where the distance between bins is zero ($i = 0$), the value $c_0$ represents the pairs residing within the exact same bin, described in Equation [](#c0-sum). 

(sec-bw-est)=
### Secondary bandwidth ($\alpha$) estimator

```{math}
:label: hat-alpha 
\hat{\alpha}(h) = 1.357 \cdot \left( \frac{\hat{S}_D(a)}{\hat{T}_D(b)} \right)^{1/7} \cdot h^{5/7}
```

Equation [](#hat-alpha) is transformed from equation [](#alpha-trans1), where both kernel function $K$ and $L$ follows Probability Density Function (PDF) of Normal Distribution. Additionally, $R(f'')$ and $R(f''')$ are estimated to $\hat{S}_D(a)$ and $\hat{T}_D(b)$ perspectively using the method of @jonesUsingNonstochasticTerms1991 with the heuristic pilot bandwidth $a$ and $b$. The practical computation of $\hat{S}_D(a)$ and $\hat{T}_D(b)$ details in [](#fist-pilot-bw) and [](#sec-pilot-bw) correspondently.

(fist-pilot-bw)=
#### The first pilot bandwidth

```{math}
:label: st-pilot-bw
a = 0.920 \cdot \hat{\lambda}N^{-1/7}
```

By Equation [](#alpha-eq), assuming both the kernel function $L$ and the true density function $f$ follow normal distributions, the heuristic bandwidth $a$ for bias optimization is obtained as Equation [](#st-pilot-bw) with $\hat{\lambda}$ is the estimated interquartile.

```{math}
:label: s-hat
\hat{S}_D(a) = \frac{1}{N(N-1)a^5} \sum_{i=0}^{N-1} \left[ c_i \cdot \phi^{\text{iv}}\left(\frac{d \cdot i}{a}\right) \right]
```

The pilot roughness estimator $\hat{S}_D(a)$ of Equation [](#hat-alpha) is estimated following Equation [](#s-hat), which originates from Equation [](#s-diagonal). Practically, due to using bin-step density ($c_i$) with the bin width $d$ calculated while dividing bins, it reduces complexity from $O(k^2)$ ($k$ is the data size) to $O(N)$ comparing to the theoretical formula.

```{math}
:label: phi-iv
\phi^{\text{iv}}(t) = \frac{1}{\sqrt{2\pi}} e^{-t^2/2} \left( t^4 - 6t^2 + 3 \right) 
```

In Equation [](#s-hat), $\phi^{\text{iv}}(t)$ is the 4-th derivative of the Normal Distribution PDF. The expansion using the [Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials) is shown by Equation [](#phi-iv). 

(sec-pilot-bw)=
#### The second pilot bandwidth

```{math}
:label: sec-pilot-bw-eq
b = 0.912 \cdot \hat{\lambda}N^{-1/9}
```

Simmilarly, $R(f''')$ is estimated and canceled out bias by using specific bandwidth $b$. Subsequently, assuming both kernel functions follow normal scale model, the heuristic bandwidth $b$ is computed following Equation [](#sec-pilot-bw-eq).  

```{math}
:label: t-hat
\begin{cases}
\begin{aligned}

\hat{T}_D(b) = \frac{1}{N(N-1)b^7} \sum_{i=0}^{N-1} \left[ c_i \cdot \phi^{\text{vi}}\left(\frac{d \cdot i}{b}\right) \right] \\
\phi^{\text{vi}}(t) = \frac{1}{\sqrt{2\pi}} e^{-t^2/2} \left(t^6 - 15t^4 + 45t^2 - 15\right)

\end{aligned}
\end{cases}
```

Likewise to [the first pilot bandwidth](#fist-pilot-bw), the component $\hat{T}_D(b)$ estimated for $R(f''')$ in Equation [](#hat-alpha) is computed by Equation [](#t-hat). Note that $\phi^{\text{vi}}(t)$ is the 6-th derivative of the PDF of Normal Distribution. 

(prim-bw-h-comp)=
### Primary bandwidth ($h$) computation

```{math}
:label: ultimate-bw
\left[  2N\sqrt{\pi} \hat{S}_D( \hat{\alpha}(h) ) \right]^{-1/5} - h = 0
```

Ultimately, the optimal bandwidth $h$ is found by Equation [](#ultimate-bw), which originates from Equation [](#amise-dev) with $K$ following normal scale model and the substitution of $R(f'')$ by $\hat{S}_D( \hat{\alpha}(h) )$. Noticeably, data size $n$ is replaced by the number of bins $N$. 

Finally, Equation [](#ultimate-bw) is solved using Brent's method [@brent2013algorithms] to find optimized bandwidth $h$. 

## Summary

Eventhough @sheatherReliableDataBasedBandwidth1991 stated that the final bandwidth $h$ is solved using Newton-Raphson method, Brent's method is employed in practice to find root of Equation [](#ultimate-bw) ensuring convergence and computational robustness.  

The biggest downside of the method is assuming real density function $f$ following another normal-scaled kernel estimator.