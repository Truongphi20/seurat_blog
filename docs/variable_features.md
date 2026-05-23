---
numbering:
    math: true
---

# Feature selection

## Introduction

## Methods

### Variance Stabilizing Transformation (`vst`)

This method computes standardized variance based on the layer input (the `count` matrix according the tutorial).

```{math}
:label: compute-var
\sigma_{i}^2 = \frac{\sum{(c_{ij} - \mu_{i})^2}}{N-1}
```

Firsly, the variance ($\sigma_{i}^2$) corresponding for each gene is computed following the formula [](#compute-var) with $c_{ij}$ is a element of gene $i$ and cell $j$ in the count matrix; $\mu_{i}$ is the mean of expression count of gene $i$; and N is the number of cell.

Subsequently, the expected variance is estimated by the Local Polynomial Regression model [@cleveland2017local], which applies linear regression to determine likelihood variance on polinomial graph in a window along the mean values. When polynomial model is identified, expected variance is returned by the mean values.         

### Mean variance plot (`mvp`)

### Dispersion (`disp`)


## Summary