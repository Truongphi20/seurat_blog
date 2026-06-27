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

```{math}
\begin{cases}
\begin{aligned}

j&: i \to 1 \\
\Delta c_j &= b_i \times b_{(i-j)}

\end{aligned}
\end{cases}
```

```{math}
c_0 = \frac{1}{2}\sum_{i=0}^{N}{b_i (b_i-1)}
```

## Summary