library(dplyr)
library(Seurat)
library(patchwork)

# commands/seurat-5.5.0/R/preprocessing5.R:66
FindVariableFeatures.StdAssay <- function(
  object,
  method = NULL,
  nfeatures = 2000L,
  layer = NULL,
  span = 0.3,
  clip = NULL,
  key = NULL,
  verbose = TRUE,
  selection.method = 'vst',
  ...
){
    layer <- "counts"
    method <- VST
    key <- 'vst'

    layer <- Layers(object = object, search = layer)
    data <- LayerData(object = object, layer = layer[1], fast = TRUE)

    hvf.info <- Seurat:::FindVariableFeatures.default(
      object = data,
      method = method,
      nfeatures = nfeatures,
      span = span,
      clip = clip,
      verbose = verbose,
      ...
    )

    colnames(x = hvf.info) <- paste(
      'vf',
      key,
      layer[1],
      colnames(x = hvf.info),
      sep = '_'
    )

    rownames(x = hvf.info) <- Features(x = object, layer = layer[1])
    object[["var.features"]] <- NULL
    object[["var.features.rank"]] <- NULL
    object[[names(x = hvf.info)]] <- NULL
    object[[names(x = hvf.info)]] <- hvf.info

    VariableFeatures(object) <- VariableFeatures(object, nfeatures=nfeatures,method = key)
    return(object)
}

# commands/seurat-5.5.0/R/preprocessing.R:4595
FindVariableFeatures.Seurat <- function(
  object,
  assay = NULL,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = 'auto',
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  ...
) {
    assay <- "RNA"
    assay.data <- FindVariableFeatures.StdAssay(
        object = object[[assay]],
        selection.method = selection.method,
        loess.span = loess.span,
        clip.max = clip.max,
        mean.function = mean.function,
        dispersion.function = dispersion.function,
        num.bin = num.bin,
        binning.method = binning.method,
        nfeatures = nfeatures,
        mean.cutoff = mean.cutoff,
        dispersion.cutoff = dispersion.cutoff,
        verbose = verbose,
        ...
    )
    object[[assay]] <- assay.data
    return(object)
}

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


## Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

## Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures.Seurat(pbmc, selection.method = "vst", nfeatures = 2000)

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