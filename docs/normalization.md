# Normalization

## Introduction

In scRNA-seq experiments, stochastic technical factors (e.g., purification, reverse transcription, and sequencing) introduce non-biological variation in cell sequencing depth. Since this noise obscures true biological signals, normalization is a critical preprocessing step to scale raw counts and enable meaningful cell-to-cell comparisons [@heumosBestPracticesSinglecell2023].

:::::{tip} Seurat command
:class: dropdown
:open: true

::::{tab-set}
:::{tab-item} Relative counts
```R
pbmc <- NormalizeData(
    pbmc, 
    normalization.method = "RC", 
    scale.factor = 10000
)
```
:::

:::{tab-item} Log normalize
```R
pbmc <- NormalizeData(
    pbmc, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
)
```
:::

:::{tab-item} Centered log ratio transformation
```R
pbmc <- NormalizeData(pbmc, normalization.method = "CLR")
```
:::

::::
:::::

## Methods

There are three methods to technically normalize data in Seurat using the `NormalizeData()` command, which were used depends on the analysis purpose.

```{image} ./static/count_matrix.png
:alt: Count matrix
:width: 100%
:align: center
```

### Relative counts (`RC`) 

$$\text{RC}_{ij} = \frac{C_{ij}}{\sum_{k=1}^{M}{C_{kj}}} \times F$$


Relative count $\text{RC}_{ij}$ represents the expression proportion of a gene $i$ within an individual cell $j$. To improve interpretability, a scale factor $F$ (defaulting to 10,000) is integrated into the calculation.    

(log-normalize)=
### Log normalize (`LogNormalize`) 

```{math}
:label: log-norm
\text{LN}_{ij} = \ln \left( \text{RC}_{ij} + 1 \right)
```

Hence the highly expressed genes possess mean, and variance much larger than low-expression genes [@ahlmann-eltzeComparisonTransformationsSinglecell2023]. To avoid data skewing in individual cell, a logarithmic transformation is employed.  

Therefore logarithm transform is applied to convert multiplicative biological relationships into an additive scale, making fold-change comparisons far more effective. Noticeably, the addition of a pseudo-count ($+1$) serves as a mathematical guardrail to prevent undefined values for zeros.

### Centered log ratio transformation (`CLR`)

$$
\begin{equation*}
\begin{cases}
    X_{ij} = C_{ij} + 1 \\[1ex]
    \text{CLR}_{ij} = \ln \left[ 1 + \dfrac{C_{ij}}{\left(\prod_{k=1}^{M} X_{kj}\right)^{\frac{1}{M}}} \right]
\end{cases}
\end{equation*}
$$

While "Relative counts" and "Log normalize" introduce the proportion of gene expression in a cell, "Centered log ratio transformation" quantifies the scaled value by the geometric mean, which acts as the center data point. 

$X_{ij}$ is the shifted count guarding zero values when calculating the geometric mean (the denominator). The $\ln(x+1)$ layer is applied to compress the dynamic range and reduce the dominance of highly expressed genes.

Note that CLR can also be applied horizontally (across cells for an individual gene) by setting the `margin` value to 2.

## Summary

| Method   |   Purpose    |   A-code    |
| :------- | :----------- | :---------- |
| Relative counts | Converts to gene expression proportion within a cell | [rc_normalise.R](https://raw.githubusercontent.com/Truongphi20/seurat_blog/refs/heads/main/algorithm_code/rc_normalise.R) |
| Log normalize (*default*)  | Compresses the dynamic range to reduce the dominance of high-expression genes in relative counts | [log_normalise.R](https://raw.githubusercontent.com/Truongphi20/seurat_blog/refs/heads/main/algorithm_code/log_normalise.R) |
| Centered log ratio transformation | Measures the scaled value by the geometric mean  | [clr_normalise.R](https://raw.githubusercontent.com/Truongphi20/seurat_blog/refs/heads/main/algorithm_code/clr_normalise.R) |


Besides three above methods which force the evenness across cells in normalised data, the SCTransform is introduced as an alternative approach by conducting statistic regression [@choudharyComparisonEvaluationStatistical2022].