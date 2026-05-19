library(dplyr)
library(Seurat)
library(patchwork)

# commands/seurat-5.5.0/src/data_manipulation.cpp:305
SparseRowVarStd_R <- function(mat, mu, sd, vmax, display_progress = FALSE) {
  if (display_progress) {
    message("Calculating feature variances of standardized and clipped values")
  }
  
  n_cells <- ncol(mat)
  n_genes <- nrow(mat)
  
  # Extract standard dgCMatrix slots
  x_vals <- mat@x
  p_ptr  <- mat@p
  
  # Pre-allocate output vector for standardized variances
  allVars <- numeric(n_genes)
  
  # Loop over each gene
  for (k in 1:n_genes) {
    # If standard deviation is 0, skip to avoid division-by-zero (variance remains 0)
    if (is.na(sd[k]) || sd[k] == 0) {
      allVars[k] <- 0
      next
    }
    
    # Get internal column pointer boundaries
    start_idx <- p_ptr[k] + 1
    end_idx   <- p_ptr[k + 1]
    
    # Catch cases where there are absolutely no non-zero elements for this feature
    if (is.na(start_idx) || is.na(end_idx) || (start_idx > end_idx)) {
      n_nonzero <- 0
    } else {
      n_nonzero <- end_idx - start_idx + 1
    }
    
    nZero <- n_cells - n_nonzero
    colSum <- 0
    
    if (n_nonzero > 0) {
      # Grab the non-zero raw counts for this gene
      gene_nonzero_vals <- x_vals[start_idx:end_idx]
      
      # Standardize the non-zero values: (value - mu) / sd
      standardized_nonzero <- (gene_nonzero_vals - mu[k]) / sd[k]
      
      # Clip values using the vmax threshold: std::min(vmax, value)
      # pmin() processes the vector element-wise against the scalar vmax
      clipped_nonzero <- pmin(vmax, standardized_nonzero)
      
      # Sum of squared deviations for non-zero items
      colSum <- sum(clipped_nonzero^2)
    }
    
    # Standardize the structural zeros: (0 - mu) / sd
    standardized_zero <- (0 - mu[k]) / sd[k]
    
    # Clip the zero-value representation as well
    clipped_zero <- pmin(vmax, standardized_zero)
    
    # Add the mathematical contribution of the omitted zeros
    colSum <- colSum + (clipped_zero^2) * nZero
    
    # Calculate sample variance of the standardized, clipped entries
    allVars[k] <- colSum / (n_cells - 1)
  }
  
  return(allVars)
}

# commands/seurat-5.5.0/src/data_manipulation.cpp:278
SparseRowVar2_R <- function(mat, mu, display_progress = FALSE) {
    
  n_cells <- ncol(mat)
  n_genes <- nrow(mat)
  
  # Extract standard dgCMatrix slots
  x_vals <- mat@x
  p_ptr  <- mat@p
  
  # Pre-allocate output vector for variances
  allVars <- numeric(n_genes)
  
  # Loop over each gene
  for (k in 1:n_genes) {
    # 0-indexed adjustment for R's 1-indexed vectors
    start_idx <- p_ptr[k] + 1
    end_idx   <- p_ptr[k + 1]
    
    # Catch cases where there are absolutely no non-zero elements for this feature
    if (is.na(start_idx) || is.na(end_idx) || (start_idx > end_idx)) {
      allVars[k] <- 0
      next
    }
    
    n_nonzero <- end_idx - start_idx + 1
    nZero <- n_cells - n_nonzero
    
    # Grab just the non-zero raw counts for this gene
    gene_nonzero_vals <- x_vals[start_idx:end_idx]
    
    # Sum of squared deviations for non-zero items: (value - mu_k)^2
    colSum <- sum((gene_nonzero_vals - mu[k])^2)
    
    # Add the mathematical contribution of the structural zeros: (0 - mu_k)^2 * nZero
    colSum <- colSum + (mu[k]^2) * nZero
    
    # Calculate sample variance (divide by N - 1)
    allVars[k] <- colSum / (n_cells - 1)
  }
  
  return(allVars)
}

# commands/seurat-5.5.0/R/preprocessing5.R:542
VST.dgCMatrix <- function(
  data,
  margin = 1L,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  verbose = TRUE,
  ...
) {
    nfeatures <- nrow(x = data)
    hvf.info <- EmptyDF(n = nfeatures)
    # Calculate feature means
    hvf.info$mean <- Matrix::rowMeans(x = data)
    # Calculate feature variance
    hvf.info$variance <- SparseRowVar2_R(
        mat = data,
        mu = hvf.info$mean,
        display_progress = FALSE
    )
    hvf.info$variance.expected <- 0L
    not.const <- hvf.info$variance > 0
    fit <- loess(
        formula = log10(x = variance) ~ log10(x = mean),
        data = hvf.info[not.const, , drop = TRUE],
        span = span
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    hvf.info$variance.standardized <- SparseRowVarStd_R(
        mat = data,
        mu = hvf.info$mean,
        sd = sqrt(x = hvf.info$variance.expected),
        vmax = clip %||% sqrt(x = ncol(x = data)),
        display_progress = verbose
    )
    # Set variable features
    hvf.info$variable <- FALSE
    hvf.info$rank <- NA
    vf <- head(
        x = order(hvf.info$variance.standardized, decreasing = TRUE),
        n = nselect
    )
    hvf.info$variable[vf] <- TRUE
    hvf.info$rank[vf] <- seq_along(along.with = vf)
    return(hvf.info)
}

# commands/seurat-5.5.0/R/preprocessing5.R:27
FindVariableFeatures.default <- function(
  object,
  method = VST,
  nfeatures = 2000L,
  verbose = TRUE,
  selection.method = selection.method,
  ...
) {
    var.gene.ouput <- VST(
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
    layer <- "counts"
    method <- VST
    key <- 'vst'

    layer <- Layers(object = object, search = layer)
    data <- LayerData(object = object, layer = layer[1], fast = TRUE)

    hvf.info <- FindVariableFeatures.default(
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