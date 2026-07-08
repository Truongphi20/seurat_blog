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

:::{tip} Why don't just solve the minimization of AMISE to find bandwidth $h$?

```{math}
:label: amise-eq
\text{AMISE}(h) = \underbrace{(nh)^{-1} R(K) \vphantom{\frac{1}{4}}}_{\text{variance}} + \underbrace{\frac{1}{4} h^4\sigma^4_K R(f'')}_{\text{bias}}
```

Asymtotic mean integrated squared error (AMISE) measures disimilarity between the beneath real model and the estimated one from the observed data. 

As mentioned by @sheatherReliableDataBasedBandwidth1991, AMISE is shown by Equation [](#amise-eq) with two components representing the variance of estimation and the bias comparing to underlying model. Particularly:

- $n$: The length of data
- $h$: The bandwidth
- $R(g)$: The roughness function ($\int{g(x)^2dx}$)
- $K$: The estimated model (Normal Distribution)
- $f$: The real underlying model

The ultimate target of the SJ method is choosing $h$ to as minimize AMISE as possible. Although the variance is straightforward for solving $h$ using mathematic methods, $R(f'')$ in the bias term obstructs the direct minimization by the unknown real model $f$, which is independent on $h$, and the integral layer.

Therefore, the whole workflow is firstly estimating $R(f'')$ and finding $h$ by minimalizing AMISE.

:::

:::{tip} How to estimate $R(f'')$?

The kernel-based method of @jonesUsingNonstochasticTerms1991 is utilized to estimate $R(f'')$.

```{math}
:label: roughness-define
R(f'') = \int_{-\infty}^{+\infty} [f''(x)]^2 dx = \int_{-\infty}^{+\infty} f^{\text{iv}}(x)f(x) dx = \mathbb{E}[f^{\text{iv}}(X)]
```

First, the roughness function of $f''$ is defined and performed integration by parts to be the expectation of 4-th derivative of $f$, which shown in Equation [](#roughness-define).  

```{math}
:label: real-f-est
\hat{f}_\alpha^{\text{iv}}(x) = \frac{1}{n\alpha^5} \sum_{j=1}^n L^{\text{iv}}\left(\frac{x - X_j}{\alpha}\right)
```

The method assumes that the real model following another [Kernel Density Estimation (KDE)](https://en.wikipedia.org/wiki/Kernel_density_estimation), which is different to the KDE when expanding AMISE (kernel function $K$, bandwidth $h$), with the kernel function $L$ and bandwidth $\alpha$. The 4-th derivative of real model estimation $\hat{f}_\alpha^{\text{iv}}(x)$ at a random point follows Equation [](#real-f-est). 


```{math}
:label: s-diagonal
\hat{S}_{\text{D}}(\alpha) = \frac{1}{n} \sum_{i=1}^n \hat{f}_\alpha^{\text{iv}}(X_i) = \frac{1}{n^2 \alpha^5} \sum_{i=1}^n \sum_{j=1}^n L^{\text{iv}}\left(\frac{X_i - X_j}{\alpha}\right)
```

To be general, the average of estimation ($\hat{S}_{\text{D}}(\alpha)$) across the data represents the validly estimated $f^{\text{iv}}(X)$, which is shown by Equation [](#s-diagonal). 

```{math}
:label: sep-s-diag
\hat{S}_{\text{D}}(\alpha) = \frac{L^{\text{iv}}(0)}{n \alpha^5} + \hat{S}_{\text{ND}}(\alpha)
```

Separating $\hat{S}_{\text{D}}(\alpha)$ into 2 parts shown by Equation [](#sep-s-diag) including the diagonal term, where $i=j$, and non-diagonal $\hat{S}_{\text{ND}}(\alpha)$, which facilitates for asymptotic expansion of expectation. 

```{math}
:label: exp-s-diag
\mathbb{E}[\hat{S}_{\text{D}}(\alpha)] \approx R(f'') + \underbrace{\frac{L^{\text{iv}}(0)}{n \alpha^5}}_{\text{Positive Bias}} - \underbrace{\frac{1}{2}\alpha^2 \sigma_L^2 R(f''')}_{\text{Negative Bias}}
```

Equation [](#exp-s-diag) applies an expectation layer to obtain $\mathbb{E}[\hat{S}_{\text{D}}(\alpha)]$, which is the ultimate estimation for $\mathbb{E}[f^{\text{iv}}(X)]$ representing the $R(f'')$ proven by Equation [](#roughness-define). After asymptotically expanding $\mathbb{E}[\hat{S}_{\text{ND}}(\alpha)]$, there are approximally two components seen as bias from the true roughness $R(f'')$.      

```{math}
:label: alpha-eq
\alpha = \left( \frac{2 L^{\text{iv}}(0)}{\sigma_L^2} \right)^{1/7} R^{-1/7}(f''') n^{-1/7}
```

To correct the estimation, the bias components is canceled out by the bandwidth $\alpha$ following Equation [](#alpha-eq), where possitive bias equals negative one. 

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

In Equation [](#t-hat), $\phi^{\text{vi}}(t)$ is the 6-th derivative of the Probability Density Function (PDF) of Normal Distribution. The expansion using the [Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials) is shown by Equation [](#six-dev-nb).  

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

### Bandwidth computation

```{math}
:label: ultimate-bw
\left[  2N\sqrt{\pi} \hat{S}_D( \hat{\alpha}(h) ) \right]^{-1/5} - h = 0
```

## Summary