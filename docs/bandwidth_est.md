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

### Bin densing

Firstly, data is divided into $N$ bins (default is 1000), and counted density.  

```{math}
:label: delta-cj

c_i = \sum_{k=i}^{N}{b_{k}b_{(k-i)}}
``` 

```{math}
:label: swapping-c0
c_1 = \sum_{k=1}^{N}{\frac{b_k (b_k-1)}{2}}
```

Separately, $c_1$ is derived by Equation [](#swapping-c0). 

## Summary