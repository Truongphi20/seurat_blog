library(dplyr)
library(Seurat)
library(patchwork)
library(pbapply)

# commands/seurat-5.5.0/R/preprocessing.R:5837
CustomNormalize <- function(data, custom_function, margin, verbose = TRUE) {
    norm.data <- apply(
        X = data,
        MARGIN = margin,
        FUN = custom_function
    )

    norm.data = Matrix::t(x = norm.data)
    colnames(x = norm.data) <- colnames(x = data)
    rownames(x = norm.data) <- rownames(x = data)
    return(norm.data)
}

# commands/seurat-5.5.0/R/preprocessing.R:4912
NormalizeData.V3Matrix <- function(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  margin = 1,
  block.size = NULL,
  verbose = TRUE,
  ...
){
    normalized.data = CustomNormalize(
        data = object,
        custom_function = function(x) {
            return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
        },
        margin = margin,
        verbose = verbose
        # across = across
    )

    return(normalized.data)
}

# commands/seurat-5.5.0/R/preprocessing5.R:311
NormalizeData.StdAssay <- function(
  object,
  normalization.method = 'LogNormalize',
  scale.factor = 1e4,
  margin = 1L,
  layer = 'counts',
  save = 'data',
  verbose = TRUE,
  ...
){
    olayer <- layer <- unique(x = "counts")
    layer <- Layers(object = object, search = "counts")
    LayerData(
      object = object,
      layer = "data",
      features = Features(x = object, layer = "counts"),
      cells = Cells(x = object, layer = "counts")
    ) <- NormalizeData.V3Matrix(
      object = LayerData(object = object, layer = "counts", fast = NA),
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      margin = margin,
      verbose = verbose,
      ...
    )
    return(object)
}

# commands/seurat-5.5.0/R/preprocessing.R:5056
NormalizeData.Seurat <- function(
  object,
  assay = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  margin = 1,
  verbose = TRUE,
  ...
) {
    assay = DefaultAssay(object = object)
    assay.data <- NormalizeData(
        object = object[[assay]],
        normalization.method = normalization.method,
        scale.factor = scale.factor,
        verbose = verbose,
        margin = margin,
        ...
    )
    object[[assay]] <- assay.data
    return(object)
}


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc <- NormalizeData.Seurat(pbmc, normalization.method = "CLR", scale.factor = 10000)

## Checking the result 
# Extract the full matrix
mat <- pbmc@assays$RNA$data

# Find rows (genes) and columns (cells) with the highest expressions
top_genes <- order(rowSums(mat > 0), decreasing = TRUE)[1:5]
top_cells <- order(colSums(mat > 0), decreasing = TRUE)[1:5]

# Subset the matrix using these top indices
dense_chunk <- mat[top_genes, top_cells]

# Print it as a standard, readable matrix
print(as.matrix(dense_chunk))

# Expected output 
#        CCAGTCTGCGGAGA-1 TTACTCGAACGTTG-1 AGAGGTCTACAGCT-1 GCGAAGGAGAGCTT-1
# TMSB4X        2.3754964        1.7581948         1.656893        1.7841300
# MALAT1        0.9396948        0.5482707         0.744437        0.5596129
# B2M           1.2728471        1.5504676         1.482067        1.3958636
# RPL13A        2.0440508        1.6537031         1.539614        1.1798928
# RPL10         1.6100127        1.7631981         1.306057        1.0520122
#        GGCACGTGTGAGAA-1
# TMSB4X        1.8699463
# MALAT1        0.7062159
# B2M           1.1156693
# RPL13A        1.6224125
# RPL10         1.6372291