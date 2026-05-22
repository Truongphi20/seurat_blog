library(dplyr)
library(Seurat)
library(patchwork)

# commands/seurat-5.5.0/src/data_manipulation.cpp:305
SparseRowVarStd_R <- function(mat, mu, sd, vmax) {
  
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
SparseRowVar2_R <- function(mat, mu) {
  
  # mat is a dgCMatrix (genes/features as rows, cells/samples as columns)
  N <- ncol(mat) # Equivalent to mat.rows() in the C++ loop after transpose
  
  # Square the non-zero elements of the sparse matrix
  mat_squared <- mat
  mat_squared@x <- mat_squared@x^2
  
  # Calculate the row sums of the squared elements: sum(x^2)
  sum_x_squared <- Matrix::rowSums(mat_squared)
  
  # Calculate the sum of the non-zero elements: sum(x)
  sum_x <- Matrix::rowSums(mat)
  
  # Apply the algebraic expansion of variance: 
  # sum((x - mu)^2) = sum(x^2) - 2 * mu * sum(x) + N * mu^2
  # This automatically accounts for implicit zeros!
  sum_squares_centered <- sum_x_squared - (2 * mu * sum_x) + (N * (mu^2))
  
  # Divide by (N - 1) for sample variance
  allVars <- sum_squares_centered / (N - 1)
  
  return(allVars)
}

# commands/seurat-5.5.0/R/preprocessing5.R:542
VST.dgCMatrix <- function(data, nselect = 2000L) {

  {nfeatures <- nrow(x = data)
  hvf.info <- EmptyDF(n = nfeatures)
  # Calculate feature means
  hvf.info$mean <- Matrix::rowMeans(x = data)
  # Calculate feature variance
  hvf.info$variance <- SparseRowVar2_R(
      mat = data,
      mu = hvf.info$mean
  )
  hvf.info$variance.expected <- 0L
  not.const <- hvf.info$variance > 0
  fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, , drop = TRUE],
      span = 0.3
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  hvf.info$variance.standardized <- SparseRowVarStd_R(
      mat = data,
      mu = hvf.info$mean,
      sd = sqrt(x = hvf.info$variance.expected),
      vmax = sqrt(x = ncol(x = data))
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
  return(hvf.info)}
}

# commands/seurat-5.5.0/R/preprocessing5.R:66
FindVariableFeatures.StdAssay <- function(object, nfeatures = 2000L){

  layer <- "counts"
  key <- 'vst'

  browser()
  data <- LayerData(object = object, layer = layer, fast = TRUE)

  hvf.info <- VST.dgCMatrix(
    data = data,
    nselect = nfeatures
  )
  rownames(x = hvf.info) <- rownames(x = data)

  colnames(x = hvf.info) <- paste(
    'vf',
    key,
    layer,
    colnames(x = hvf.info),
    sep = '_'
  )

  rownames(x = hvf.info) <- Features(x = object, layer = layer)
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