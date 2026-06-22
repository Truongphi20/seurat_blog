# SCTransform workflow

## Introduction

:::{tip} Seurat command

```R
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
```

:::

## Workflow

### Fitting model

At the begining, a number of genes is randomly selected (default is 2000) depending on quantity of expression. They are divided into batches containing 500 genes and fitting Gamma-Poisson general linear model (Gamma-Poisson GLM) separately. See [](./glmGamPoi.md) for more detail about fitting model.   

### Regularizing model

### Pearson residuals

## Summary