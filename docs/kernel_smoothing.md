---
title: Kernel Regression Smoother
numbering:
    math: true
---

## Introduction

:::{tip} Seurat command
This command is captured from the process of [SCTransform](./sctransform.md).

```R
# https://github.com/satijalab/sctransform/blob/v0.4.3/R/vst.R#L857-L858
model_pars_fit[o, 'dispersion_par'] <- 
    ksmooth(
        x = genes_log_gmean_step1, 
        y = model_pars[, 'dispersion_par'],
        x.points = x_points, 
        bandwidth = bw, 
        kernel='normal'
    )$y
```

Command explanation:

- `x`: an attributes of data in increasing order.
- `y`: fitted values corresponding to `x`.
- `x.points`: where evaluating smoothed fit.
- `kernel`: the kernel function is ultilized.
- `bandwidth`: the scaling factor used for kernel function. 

:::

## Workflow

## Summary