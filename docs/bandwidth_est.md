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

### Bin-width density

Firstly, data is divided into $N$ bins (default is 1000), and counted density for each bin ($b_k, 0 \le k \le N-1$), which is prepared for identifying total denstity for each bin-width distance.  

```{math}
:label: ci-sum

c_i = \sum_{k=i}^{N-1}{b_{k}b_{(k-i)}}
``` 

Equation [](#ci-sum) measures total density of unique pairs of points ($c_i$) between 2 bins being far from each other $i$ bin(s). 

```{math}
:label: c0-sum
c_0 = \sum_{k=0}^{N-1}{ b_k \choose 2 } = \sum_{k=0}^{N-1}{\frac{b_k (b_k-1)}{2}}
```

When the distance between bins (bin-width) is 0, $c_0$ is computed by selecting the unique pairs of points from the same bins, described in Equation [](#c0-sum). 

## Summary