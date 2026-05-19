library(dplyr)
library(Seurat)
library(patchwork)

# commands/seurat-5.5.0/R/preprocessing5.R:704
DISP <- function(
  data,
  nselect = 2000L,
  verbose = TRUE,
  ...
) {
  hvf.info <- CalcDispersion(object = data, verbose = verbose, ...)
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA
  vf <- head(
    x = order(hvf.info$mvp.dispersion, decreasing = TRUE),
    n = nselect
  )
  hvf.info$variable[vf] <- TRUE
  hvf.info$rank[vf] <- seq_along(along.with = vf)
  return(hvf.info)
}

# commands/seurat-5.5.0/R/preprocessing5.R:1832
MVP <- function(
  data,
  verbose = TRUE,
  nselect = 2000L,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  ...
){
    hvf.info <- DISP(data = data, nselect = nselect, verbose = verbose)
    hvf.info$variable <- FALSE
    hvf.info$rank <- NA
    hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), , drop = FALSE]
    means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
    dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
    hvf.info[which(x = means.use & dispersions.use), 'variable'] <- TRUE
    rank.rows <- rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    selected.indices <- which(rownames(x = hvf.info) %in% rank.rows)
    hvf.info$rank[selected.indices] <- seq_along(selected.indices)
    hvf.info <- hvf.info[order(as.numeric(row.names(hvf.info))), ]

    return(hvf.info)
}

# commands/seurat-5.5.0/R/preprocessing5.R:27
FindVariableFeatures.default <- function(
  object,
  nfeatures = 2000L,
  verbose = TRUE,
  selection.method = selection.method,
  ...
){
    var.gene.ouput <- Seurat:::MVP(
        data = object,
        nselect = nfeatures,
        verbose = verbose,
        ...
    )
    rownames(x = var.gene.ouput) <- rownames(x = object)
    return(var.gene.ouput)
}

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
    layer <- "data"
    key <- 'mvp'

    layer <- Layers(object = object, search = layer)
    data <- LayerData(object = object, layer = layer[1], fast = TRUE)

    hvf.info <- FindVariableFeatures.default(
      object = data,
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
pbmc <- FindVariableFeatures.Seurat(pbmc, selection.method = "mean.var.plot", nfeatures = 2000)

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