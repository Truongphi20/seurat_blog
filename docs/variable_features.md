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

Firsly, the squared variance ($\sigma_{i}^2$) corresponding for each gene is computed following the formular [](#compute-var) with $c_{ij}$ is a element of gene $i$ and cell $j$ in the count matrix; $\mu_{i}$ is the mean of expression count of gene $i$; and N is the number of cell.

### Mean variance plot (`mvp`)

### Dispersion (`disp`)


## Summary