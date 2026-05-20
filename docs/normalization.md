# Normalization

## Introduction

In scRNA-seq experiments, stochastic technical factors (e.g., purification, reverse transcription, and sequencing) introduce non-biological variation in cell sequencing depth. Since this noise obscures true biological signals, normalization is a critical preprocessing step to scale raw counts and enable meaningful cell-to-cell comparisons [@heumosBestPracticesSinglecell2023].

## Methods

There are three methods to normalise data in Seurat, which were used depends on data input.

| Name                               | Method                         |  Meaning     |
| :--------------------------------- | :----------------------------- | :----------- |
| Relative counts (`RC`)                   |    |    | 
| Log normalize   (`LogNormalize`)                   |    |    |
| Centered log ratio transformation  (`CLR`) |    |    |

## Summary