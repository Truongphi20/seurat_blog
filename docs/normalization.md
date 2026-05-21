# Normalization

## Introduction

In scRNA-seq experiments, stochastic technical factors (e.g., purification, reverse transcription, and sequencing) introduce non-biological variation in cell sequencing depth. Since this noise obscures true biological signals, normalization is a critical preprocessing step to scale raw counts and enable meaningful cell-to-cell comparisons [@heumosBestPracticesSinglecell2023].

## Methods

There are three methods to normalise data in Seurat, which were used depends on data input.

```{image} ./static/count_matrix.png
:alt: Count matrix
:width: 60%
:align: center
```

### Relative counts (`RC`) 

$$\text{RC}_{ij} = \frac{C_{ij}}{\sum_{k=1}^{M}{C_{kj}}} \times \text{Scale Factor}$$


### Log normalize (`LogNormalize`) 

$$\text{LN}_{ij} = \ln \left( \text{RC}_{ij} + 1 \right)$$

### Centered log ratio transformation (`CLR`)

$$
\begin{equation*}
\begin{cases}
    X_{ij} = C_{ij} + 1 \\[1ex]
    \text{CLR}_{ij} = \ln \left[ 1 + \dfrac{X_{ij}}{\left(\prod_{k=1}^{M} X_{kj}\right)^{\frac{1}{M}}} \right]
\end{cases}
\end{equation*}
$$

## Summary