library(dplyr)
library(Seurat)
library(patchwork)


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


## Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


## Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "mean.var.plot", nfeatures = 2000)

## Check output
print(pbmc@assays$RNA@meta.data[1:10,1:3])

## Expected output
#    vf_mvp_data_mvp.mean vf_mvp_data_mvp.dispersion
# 1           0.013246433                   1.432911
# 2           0.004588393                   1.458647
# 3           0.005542619                   1.325485
# 4           0.002583538                   0.859281
# 5           0.026815606                   1.457617
# 6           0.373210434                   1.879582
# 7           0.015307744                   1.775412
# 8           0.012649687                   1.710906
# 9           0.003562567                   1.270087
# 10          0.217553664                   1.694766
#    vf_mvp_data_mvp.dispersion.scaled
# 1                         -0.6584199
# 2                         -0.6080841
# 3                         -0.8685312
# 4                         -1.7803661
# 5                         -0.6100982
# 6                         -0.4387589
# 7                          0.0114675
# 8                         -0.1146983
# 9                         -0.9768836
# 10                        -0.1462672