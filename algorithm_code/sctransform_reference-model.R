library(Seurat)
library(ggplot2)
library(sctransform)
library(Azimuth)


pbmc_data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc_data)

# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# Download the reference
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc", seconds = 30L)


# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt")