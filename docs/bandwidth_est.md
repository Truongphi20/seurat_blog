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

### Second pilot bandwidth

```{math}
b = 0.912 \cdot \hat{\lambda}N^{-1/9}
```

```{math}
\hat{T}_D(b) = \frac{1}{N(N-1)b^7} \sum_{i=1}^{N} \sum_{j=1}^{N} \phi^{\text{vi}}\left(\frac{X_i - X_j}{b}\right) \quad \text{}
```

```{math}
\begin{cases}
\begin{aligned}

t_{ij} &= \frac{X_i - X_j}{b} \\
\phi^{\text{vi}}(t_{ij}) &= \frac{1}{\sqrt{2\pi}} e^{-t_{ij}^2/2} \left(t_{ij}^6 - 15t_{ij}^4 + 45t_{ij}^2 - 15\right)

\end{aligned}
\end{cases}
```

## Summary