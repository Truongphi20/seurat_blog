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
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

## Check output
print(pbmc@assays$RNA@meta.data[1:10,1:3])

## Expected output
#    vf_vst_counts_mean vf_vst_counts_variance vf_vst_counts_variance.expected
# 1         0.003333333            0.003323453                     0.003575582
# 2         0.001111111            0.001110288                     0.001112798
# 3         0.001851852            0.001849107                     0.001921811
# 4         0.001111111            0.001110288                     0.001112798
# 5         0.006666667            0.006624676                     0.007342308
# 6         0.106666667            0.158310485                     0.203482316
# 7         0.003333333            0.003323453                     0.003575582
# 8         0.002592593            0.002586829                     0.002744432
# 9         0.001111111            0.001110288                     0.001112798
# 10        0.078888889            0.145311844                     0.138583325