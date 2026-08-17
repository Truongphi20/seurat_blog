---
title: "Seurat Under the Hood"
description: "Deep dives into Seurat internals"
---

## Motivation

[Seurat](https://satijalab.org/seurat/) is a popular standard bioinformatic tool to explore and analyze single cell data. The vignette (tutorial) ["Get started with Seurat"](https://satijalab.org/seurat/articles/get_started_v5_new) is seemed to be a foundation material for newcomers to get used to analysis steps in this domain.    

Eventhough the vignette provides a solid scaffold for a primaritive scRNA-seq analysis by displaying the command lines and the brief abstractions, it is still lacked of the thoughtful insight and the midset of using algorithmic and statistic methods. It is an obstacle for someone being curious about the "under-the-hood" but insufficient in technical foundation.    

I write this blog as a extension of the Seurat vignette, exploring the proceeses run beneath the "easy-going" commands and delivering their reasoning.  

## Content