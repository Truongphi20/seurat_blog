# Scaling data

## Introduction

## Methods

```{math}
:label: scaling-func
\Large
\begin{cases}
\begin{aligned}

\mu_i &= \frac{1}{N}\sum_{k=1}^{N}{RC_{ik}} \\
\sigma_i &= \sqrt{\frac{\sum_{k=1}^{N}{ (RC_{ik} - \mu_i)^2 }}{N-1}} \\
SR_{ij} &= min \left(\frac{RC_{ik} - \mu_i}{\sigma_i}, SR_{\text{max}} \right) 

\end{aligned}
\end{cases}
```

## Summary