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

```{image} ./static/bin-dens.png
```

Firstly, data is divided into $N$ bins (default is 1000), and counted density. Go through each bin ($b_i$), the former and current bins (predecessor bins) are weighted reordered with the weight being current density.   

```{math}
:label: delta-cj

\begin{cases}
\begin{aligned}

j&: i \to 1 \\
\Delta c_j &= b_i \times b_{(i-j)}

\end{aligned}
\end{cases}
```

Particularly, the changing $\Delta c_j$ of new value $c_j$ of each predecessor bin following Equation [](#delta-cj). At the $b_i$ stop, the weighted swappings ($b_0$ is included) occur from $b_i$ to $b_1$.   

```{math}
:label: swapping-c0
c_0 = \sum_{i=0}^{N}{\frac{b_i (b_i-1)}{2}}
```

Separately, $c_0$ is derived by Equation [](#swapping-c0). 

## Summary