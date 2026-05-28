# Scaling data

## Introduction

## Methods

Scaled-row value ($SR_{ij}$) of each gene $i$ and cell $j$, which is scaled vertically by gene across cells, being from the normalized matrix ([log-normalize matrix](#log-normalize) in the case of tutorial).

```{math}
:label: scaling-func
\Large
\begin{cases}
\begin{aligned}

\mu_i &= \frac{1}{N}\sum_{k=1}^{N}{LN_{ik}} \\
\sigma_i &= \sqrt{\frac{\sum_{k=1}^{N}{ (LN_{ik} - \mu_i)^2 }}{N-1}} \\
SR_{ij} &= min \left(\frac{LN_{ik} - \mu_i}{\sigma_i}, SR_{\text{max}} \right) 

\end{aligned}
\end{cases}
```

At the beginning, mean value of each gene ($\mu_i$) is computed by the average of log-normalized values across cells for gene $i$, with $N$ is the total number of cells. From that, standard deviation ($\sigma_i$) is calculated for each gene. The normalized values are standardized by mean and standard deviation of each gene with the scaled-max (default is 10, set by argument `scale.max`)        

## Summary