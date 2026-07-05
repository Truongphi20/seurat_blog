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

```{math}
:label: hat-alpha 
\hat{\alpha}_2(h) = 1.357 \cdot \left( \frac{\hat{S}_D(a)}{\hat{T}_D(b)} \right)^{1/7} \cdot h^{5/7}
```

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

### The second pilot bandwidth

```{math}
:label: sec-pilot-bw
b = 0.912 \cdot \hat{\lambda}N^{-1/9}
```

At the beginning, the initial heuristic bandwidth ($b$) is computed following Equation [](#sec-pilot-bw) with $\hat{\lambda}$ is the estimated interquartile. 

```{math}
:label: t-hat
\hat{T}_D(b) = \frac{1}{N(N-1)b^7} \sum_{i=0}^{N-1} \left[ c_i \cdot \phi^{\text{vi}}\left(\frac{d \cdot i}{b}\right) \right]
```

The dominator component of Equation [](#hat-alpha) is estimated following Equation [](#t-hat). Practically, due to using bin-step density ($c_i$) with the bin width $d$ calculated while dividing bins, it reduces complexity from $O(n^2)$ to $O(n)$ comparing to the theoretical formula mentioned by @sheatherReliableDataBasedBandwidth1991.   

```{math}
:label: six-dev-nb
\phi^{\text{vi}}(t) = \frac{1}{\sqrt{2\pi}} e^{-t^2/2} \left(t^6 - 15t^4 + 45t^2 - 15\right)
```

In Equation [](#t-hat), $\phi^{\text{vi}}(t)$ is the 6-th derivative of the Probability Density Function (PDF) of Normal Distribution. The expansion using the [Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials) is shown by Equation [](#six-dev-nb).  

### The first pilot bandwidth

```{math}
:label: st-pilot-bw

a = 0.920 \cdot \hat{\lambda}N^{-1/7}

```

```{math}
:label: s-hat
\hat{S}_D(a) = \frac{1}{N(N-1)a^5} \sum_{i=0}^{N-1} \left[ c_i \cdot \phi^{\text{iv}}\left(\frac{d \cdot i}{a}\right) \right]
```

```{math}
:label: four-dev-nb
\phi^{\text{iv}}(t) = \frac{1}{\sqrt{2\pi}} e^{-t^2/2} \left( t^4 - 6t^2 + 3 \right) 
```

## Summary