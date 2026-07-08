---
numbering:
    math: true
---

# Bandwidth Estimation

## Introduction

:::{tip} Seurat command
This command is captured from the process of [SCTransform](./sctransform.md).

```R
bw.SJ(x)
```

Here `x` is an arbitrary numeric vector.
:::

## Workflow

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

Summarily, $R(f'') \approx \hat{S}_{\text{D}}(\alpha(h_{\text{opt}}))$, which is performed as Equation [](s-diagonal) with an independent KDE obtaining $\alpha$ depended on the optimized bandwidth $h_{\text{opt}}$ shown as Equation [](#alpha-trans1). From that, $h_{\text{opt}}$ can be solved by the equation of subsituting $\hat{S}_{\text{D}}(\alpha(h_{\text{opt}}))$ to Equation [](#amise-dev). 

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

### Bandwidth for bias optimization 

```{math}
:label: hat-alpha 
\hat{\alpha}(h) = 1.357 \cdot \left( \frac{\hat{S}_D(a)}{\hat{T}_D(b)} \right)^{1/7} \cdot h^{5/7}
```

From Equation [](#alpha-trans1), assuming both kernel function $K$ and $L$ follows Probability Density Function (PDF) of Normal Distribution. Additionally, $R(f'')$ and $R(f''')$ are estimated to $\hat{S}_D(a)$ and $\hat{T}_D(b)$ perspectively using the method of @jonesUsingNonstochasticTerms1991 with the heuristic pilot bandwidth $a$ and $b$. The practical computation of $\hat{S}_D(a)$ and $\hat{T}_D(b)$ details in [](#fist-pilot-bw) and [](#sec-pilot-bw) correspondently.

(fist-pilot-bw)=
#### The first pilot bandwidth

```{math}
:label: st-pilot-bw
a = 0.920 \cdot \hat{\lambda}N^{-1/7}
```

By assuming both the pilot kernel and the unknown density to be normal distributions in the method of @jonesUsingNonstochasticTerms1991, the heuristic bandwidth $a$ for bias optimization is obtained following Equation [](#st-pilot-bw).

```{math}
:label: s-hat
\begin{cases}
\begin{aligned}

\hat{S}_D(a) &= \frac{1}{N(N-1)a^5} \sum_{i=0}^{N-1} \left[ c_i \cdot \phi^{\text{iv}}\left(\frac{d \cdot i}{a}\right) \right] \\
\phi^{\text{iv}}(t) &= \frac{1}{\sqrt{2\pi}} e^{-t^2/2} \left( t^4 - 6t^2 + 3 \right) 

\end{aligned}
\end{cases}
```

Likewise to [the second pilot bandwidth](#sec-pilot-bw), the numerator component $\hat{S}_D(a)$ for Equation [](#hat-alpha) is computed by Equation [](#s-hat). Note that $\phi^{\text{iv}}(t)$ is the 4-th derivative of the PDF of Normal Distribution. 

(sec-pilot-bw)=
#### The second pilot bandwidth

```{math}
:label: sec-pilot-bw-eq
b = 0.912 \cdot \hat{\lambda}N^{-1/9}
```

At the beginning, the initial heuristic bandwidth ($b$) is computed following Equation [](#sec-pilot-bw-eq) with $\hat{\lambda}$ is the estimated interquartile. 

```{math}
:label: t-hat
\hat{T}_D(b) = \frac{1}{N(N-1)b^7} \sum_{i=0}^{N-1} \left[ c_i \cdot \phi^{\text{vi}}\left(\frac{d \cdot i}{b}\right) \right]
```

The dominator component of Equation [](#hat-alpha) is estimated following Equation [](#t-hat). Practically, due to using bin-step density ($c_i$) with the bin width $d$ calculated while dividing bins, it reduces complexity from $O(k^2)$ ($k$ is the number of data points) to $O(N)$ comparing to the theoretical formula mentioned by @sheatherReliableDataBasedBandwidth1991.   

```{math}
:label: six-dev-nb
\phi^{\text{vi}}(t) = \frac{1}{\sqrt{2\pi}} e^{-t^2/2} \left(t^6 - 15t^4 + 45t^2 - 15\right)
```

In Equation [](#t-hat), $\phi^{\text{vi}}(t)$ is the 6-th derivative of the Normal Distribution PDF. The expansion using the [Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials) is shown by Equation [](#six-dev-nb).  

### Bandwidth computation

```{math}
:label: ultimate-bw
\left[  2N\sqrt{\pi} \hat{S}_D( \hat{\alpha}(h) ) \right]^{-1/5} - h = 0
```

## Summary