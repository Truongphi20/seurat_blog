library(dplyr)
library(Seurat)
library(patchwork)

# commands/seurat-5.5.0/src/data_manipulation.cpp:335
FastLogVMR <- function(mat) {
  
  ncols <- ncol(mat)
  
  # 1. Transform the non-zero values to exp(x) - 1
  mat_expm1 <- mat
  mat_expm1@x <- expm1(mat_expm1@x)
  
  # 2. Compute the means of the untransformed data
  rm_means <- Matrix::rowSums(mat_expm1) / ncols
  
  # 3. Compute the sample variances of the untransformed data
  # Mathematically, Var(X) = (Sum(X^2) - n * Mean^2) / (n - 1)
  # We square the non-zero elements directly to preserve structural sparsity
  mat_expm1_sq <- mat_expm1
  mat_expm1_sq@x <- (mat_expm1_sq@x)^2
  
  sum_squares <- Matrix::rowSums(mat_expm1_sq)
  v_variances <- (sum_squares - ncols * (rm_means^2)) / (ncols - 1)
  
  # 4. Calculate the Log VMR
  rowdisp <- log(v_variances / rm_means)
  
  return(rowdisp)
}

# commands/seurat-5.5.0/src/data_manipulation.cpp:255
FastExpMean <- function(mat) {
    
  ncols <- ncol(mat)
  
  # 1. Transform the non-zero values to exp(x) - 1
  # Operating directly on the x slot of a dgCMatrix preserves structural sparsity
  mat_expm1 <- mat
  mat_expm1@x <- expm1(mat_expm1@x)
  
  # 2. Sum up rows, divide by number of columns, and log1p the result
  row_sums <- Matrix::rowSums(mat_expm1)
  rm_means <- row_sums / ncols
  rowmeans <- log1p(rm_means)
  
  return(rowmeans)
}

# commands/seurat-5.5.0/R/preprocessing5.R:629
CalcDispersion <- function(object){

  feature.mean <- FastExpMean(object)
  feature.dispersion <- FastLogVMR(object)

  names(x = feature.mean) <- names(
  x = feature.dispersion) <- rownames(x = object)
  feature.dispersion[is.na(x = feature.dispersion)] <- 0
  feature.mean[is.na(x = feature.mean)] <- 0

  data.x.breaks <- 20
  data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks,
                  include.lowest = TRUE)
  
  names(x = data.x.bin) <- names(x = feature.mean)
  mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
  sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
  feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
      sd.y[as.numeric(x = data.x.bin)]
  names(x = feature.dispersion.scaled) <- names(x = feature.mean)
  hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)

  rownames(x = hvf.info) <- rownames(x = object)
  colnames(x = hvf.info) <- paste0('mvp.', c('mean', 'dispersion', 'dispersion.scaled'))
  
  return(hvf.info)
}

# commands/seurat-5.5.0/R/preprocessing5.R:704
DISP <- function(
  data,
  nselect = 2000L
) {
  hvf.info <- CalcDispersion(object = data)
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
MVP <- function(data, nselect = 2000L){

    mean.cutoff = c(0.1, 8)
    dispersion.cutoff = c(1, Inf)

    hvf.info <- DISP(data = data, nselect = nselect)
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
  nfeatures = 2000L
){
    var.gene.ouput <- MVP(data = object, nselect = nfeatures)
    rownames(x = var.gene.ouput) <- rownames(x = object)
    return(var.gene.ouput)
}

# commands/seurat-5.5.0/R/preprocessing5.R:66
FindVariableFeatures.StdAssay <- function(object, nfeatures = 2000L){

    layer <- "data"
    key <- 'mvp'

    layer <- Layers(object = object, search = layer)
    data <- LayerData(object = object, layer = layer[1], fast = TRUE)

    hvf.info <- FindVariableFeatures.default(object = data, nfeatures = nfeatures)

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
FindVariableFeatures.Seurat <- function(object, nfeatures = 2000) {
  assay <- "RNA"
  assay.data <- FindVariableFeatures.StdAssay(
      object = object[[assay]],
      nfeatures = nfeatures
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
pbmc <- FindVariableFeatures.Seurat(pbmc, nfeatures = 2000)

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