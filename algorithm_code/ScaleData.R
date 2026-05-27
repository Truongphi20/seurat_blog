library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

# commands/seurat-5.5.0/src/data_manipulation.cpp:148
FastSparseRowScale <- function(mat, scale = TRUE, center = TRUE, scale_max = 10) {

  scale = TRUE
  center = TRUE
  scale_max = 10
   
  # Get dimensions of original matrix
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  
  # Pre-allocate the dense matrix for the output
  scaled_mat <- matrix(0, nrow = n_rows, ncol = n_cols)
  
  # Loop through each row
  for (k in 1:n_rows) {
    
    # Extract the row as a sparse vector
    row_vals <- mat[k, , drop = TRUE] # numeric vector, mostly zeros
    
    # Count non-zero elements
    nz_indices <- which(row_vals != 0)
    nnZero <- length(nz_indices)
    
    # Calculate Mean (colMean in C++)
    # C++ loops through non-zero values to sum them, then divides by total rows of transposed (cols of original)
    row_sum <- sum(row_vals[nz_indices])
    colMean <- row_sum / n_cols
    
    # Calculate Standard Deviation
    if (scale) {
      colSdev <- 0
      if (center) {
        # Sum squared differences for non-zero elements
        if (nnZero > 0) {
          colSdev <- sum((row_vals[nz_indices] - colMean)^2)
        }
        # Add the contribution of the implicit zeros: (0 - colMean)^2 * (total_cells - nnZero)
        colSdev <- colSdev + (colMean^2) * (n_cols - nnZero)
      } else {
        # If not centering, just sum the squares of non-zero elements
        if (nnZero > 0) {
          colSdev <- sum(row_vals[nz_indices]^2)
        }
      }
      # Calculate sample standard deviation (denominator is N - 1)
      colSdev <- sqrt(colSdev / (n_cols - 1))
    } else {
      colSdev <- 1
    }
    
    # If center is FALSE, the mean subtraction is skipped (colMean set to 0)
    if (!center) {
      colMean <- 0
    }
    
    # Scale, Center, and Clip (scale_max)
    if (colSdev == 0) {
      scaled_row <- rep(0, n_cols)
    } else {
      scaled_row <- (row_vals - colMean) / colSdev
    }
    
    # Apply max clipping threshold (scale_max)
    scaled_row[scaled_row > scale_max] <- scale_max
    
    # Store into output matrix
    scaled_mat[k, ] <- scaled_row
  }
  
  return(scaled_mat)
}

# commands/seurat-5.5.0/R/preprocessing.R:5111
ScaleData.default <- function(object, features = NULL)
{
    block.size = 1000

    features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
    object <- object[features, , drop = FALSE]
    object.names <- dimnames(x = object)

    split.by <- TRUE
    split.cells <- split(x = colnames(x = object), f = split.by)

    scale.function <- FastSparseRowScale

    scaled.data <- matrix(
      data = NA_real_,
      nrow = nrow(x = object),
      ncol = ncol(x = object),
      dimnames = object.names
    )
    max.block <- ceiling(x = length(x = features) / block.size)

    for (i in 1:max.block) {
        my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
        my.inds <- my.inds[my.inds <= length(x = features)]
        arg.list <- list(mat = object[features[my.inds], split.cells[["TRUE"]], drop = FALSE])
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[my.inds], split.cells[["TRUE"]]])
        scaled.data[features[my.inds], split.cells[["TRUE"]]] <- data.scale
    }

    dimnames(x = scaled.data) <- object.names
    scaled.data[is.na(x = scaled.data)] <- 0

    return(scaled.data)
}

# commands/seurat-5.5.0/R/preprocessing5.R:359
ScaleData.StdAssay <- function(object, features = NULL)
{
    layer <- 'data'

    ldata <- LayerData(object = object, layer = layer, features = features)

    ldata <- ScaleData.default(
      object = ldata,
      features = features
    )

    LayerData(object = object, layer = 'scale.data', features = rownames(ldata)) <- ldata

    return(object)
}

# commands/seurat-5.5.0/R/preprocessing.R:5520
ScaleData.Seurat <- function(
  object,
  features = NULL
) {
    assay <- "RNA"

    assay.data <- ScaleData.StdAssay(
        object = object[[assay]],
        features = features
    )
    object[[assay]] <- assay.data
    return(object)
}


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Normalization
pbmc <- NormalizeData(pbmc)

# Find highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData.Seurat(pbmc, features = all.genes)

## Checking the result 
print(pbmc[["RNA"]]$scale.data[1:5,1:5])

# Expected output
#               AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1
# AL627309.1         -0.05744997      -0.05744997      -0.05744997
# AP006222.2         -0.03318769      -0.03318769      -0.03318769
# RP11-206L10.2      -0.04118635      -0.04118635      -0.04118635
# RP11-206L10.9      -0.03325679      -0.03325679      -0.03325679
# LINC00115          -0.08128413      -0.08128413      -0.08128413
#               AAACCGTGCTTCCG-1 AAACCGTGTATGCG-1
# AL627309.1         -0.05744997      -0.05744997
# AP006222.2         -0.03318769      -0.03318769
# RP11-206L10.2      -0.04118635      -0.04118635
# RP11-206L10.9      -0.03325679      -0.03325679
# LINC00115          -0.08128413      -0.08128413