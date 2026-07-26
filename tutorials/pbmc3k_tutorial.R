## Setup the Seurat Object
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


## Standard pre-processing workflow

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

## Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

## Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

## Determine the ‘dimensionality’ of the dataset
ElbowPlot(pbmc)

# ## Cluster the cells
# pbmc <- FindNeighbors(pbmc, dims = 1:10)
# pbmc <- FindClusters(pbmc, resolution = 0.5)

# ## Run non-linear dimensional reduction (UMAP/tSNE)
# pbmc <- RunUMAP(pbmc, dims = 1:10)

# ## Finding differentially expressed features (cluster biomarkers)

# # find all markers of cluster 2
# cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
# head(cluster2.markers, n = 5)

# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
# head(cluster5.markers, n = 5)

# # find markers for every cluster compared to all remaining cells, report only the positive
# # ones
# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
# pbmc.markers %>%
#     group_by(cluster) %>%
#     dplyr::filter(avg_log2FC > 1)
# cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# ## Assigning cell type identity to clusters
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#     "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)