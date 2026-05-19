library(dplyr)
library(Seurat)
library(patchwork)

# commands/seurat-5.5.0/src/data_manipulation.cpp:335
FastLogVMR <- function(mat, display_progress = FALSE) {
    
  n_cells <- ncol(mat)
  n_genes <- nrow(mat)
  
  # Extract standard dgCMatrix structural slots
  x_vals <- mat@x
  p_ptr  <- mat@p
  
  # Pre-allocate output vector for row dispersions
  rowdisp <- numeric(n_genes)
  
  # Loop over each gene
  for (k in 1:n_genes) {
    # Get internal column pointer boundaries for gene k
    start_idx <- p_ptr[k] + 1
    end_idx   <- p_ptr[k + 1]
    
    # Catch cases where there are absolutely no non-zero elements for this feature
    if (is.na(start_idx) || is.na(end_idx) || (start_idx > end_idx)) {
      nnZero <- 0
    } else {
      nnZero <- end_idx - start_idx + 1
    }
    
    rm_sum <- 0
    
    # --- Step 1: Calculate the mean in linear space ---
    if (nnZero > 0) {
      gene_nonzero_vals <- x_vals[start_idx:end_idx]
      rm_sum <- sum(expm1(gene_nonzero_vals))
    }
    
    # Structural zeros evaluate to 0 in expm1 space, so they don't add to the sum
    rm <- rm_sum / n_cells
    
    # Guard against completely unexpressed genes to avoid division-by-zero
    if (rm == 0) {
      rowdisp[k] <- NA
      next
    }
    
    # --- Step 2: Calculate the variance in linear space ---
    v_sum <- 0
    if (nnZero > 0) {
      # Square of deviations for non-zero entries in linear space: (expm1(x) - mean)^2
      v_sum <- sum((expm1(gene_nonzero_vals) - rm)^2)
    }
    
    # Add the contribution of the structural zero entries: (0 - mean)^2 * nZero
    nZero <- n_cells - nnZero
    v_sum <- v_sum + (nZero * (rm^2))
    
    # Calculate sample variance
    v <- v_sum / (n_cells - 1)
    
    # --- Step 3: Compute the final log(Variance-to-Mean Ratio) ---
    # Guard against 0 variance or negative numbers before taking the natural log
    if (v <= 0) {
      rowdisp[k] <- NA
    } else {
      rowdisp[k] <- log(v / rm)
    }
  }
  
  return(rowdisp)
}

# commands/seurat-5.5.0/src/data_manipulation.cpp:255
FastExpMean <- function(mat, display_progress = FALSE) {
   
  n_cells <- ncol(mat)
  n_genes <- nrow(mat)
  
  # Extract standard dgCMatrix structural slots
  x_vals <- mat@x
  p_ptr  <- mat@p
  
  # Pre-allocate output vector for row means
  rowmeans <- numeric(n_genes)
  
  # Loop over each gene (analogous to the C++ outer loop after transposing)
  for (k in 1:n_genes) {
    # Get internal column pointer boundaries for gene k
    start_idx <- p_ptr[k] + 1
    end_idx   <- p_ptr[k + 1]
    
    # Catch cases where there are absolutely no non-zero elements for this feature
    if (is.na(start_idx) || is.na(end_idx) || (start_idx > end_idx)) {
      n_nonzero <- 0
    } else {
      n_nonzero <- end_idx - start_idx + 1
    }
    
    rm_sum <- 0
    
    if (n_nonzero > 0) {
      # Grab just the non-zero counts for this gene
      gene_nonzero_vals <- x_vals[start_idx:end_idx]
      
      # Sum up the un-logged values: exp(x) - 1
      # expm1() is vector-optimized in R and matches C++ expm1()
      rm_sum <- sum(expm1(gene_nonzero_vals))
    }
    
    # The contribution of structural zeros to the sum is 0, since exp(0) - 1 = 0.
    # Therefore, we can just divide directly by total cell count (ncols in C++)
    rm_avg <- rm_sum / n_cells
    
    # Log-transform the averaged result back: log(1 + mean)
    rowmeans[k] <- log1p(rm_avg)
  }
  
  return(rowmeans)
}

# commands/seurat-5.5.0/R/preprocessing5.R:629
CalcDispersion <- function(
  object,
  mean.function = FastExpMean,
  dispersion.function = FastLogVMR,
  num.bin = 20,
  binning.method = "equal_width",
  verbose = TRUE,
  ...
){
    feature.mean <- FastExpMean(object, verbose)
    feature.dispersion <- FastLogVMR(object, verbose)

    names(x = feature.mean) <- names(
    x = feature.dispersion) <- rownames(x = object)
    feature.dispersion[is.na(x = feature.dispersion)] <- 0
    feature.mean[is.na(x = feature.mean)] <- 0

    data.x.breaks <- num.bin
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