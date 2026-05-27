# Scaling data

## Introduction

## Methods

```{math}
\Large \mu_i = \frac{1}{N}\sum_{k=1}^{N}{RC_{ik}}
```

```{math}
\Large \sigma_i = \sqrt{\frac{\sum_{k=1}^{N}{ (RC_{ik} - \mu_i)^2 }}{N-1}} 
```

```{math}
\Large SR_{ij} = min \left(\frac{RC_{ik} - \mu_i}{\sigma_i}, SR_{\text{max}} \right)
```

## Summary