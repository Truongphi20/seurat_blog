library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

# commands/seurat-5.5.0/src/data_manipulation.cpp:115
LogNorm <- function(data, scale_factor, display_progress = FALSE) {
    col_sums <- Matrix::colSums(data)
    cell_indices <- rep.int(1:ncol(data), diff(data@p))
    data@x <- log1p(data@x / col_sums[cell_indices] * scale_factor)
    return(data)
}

# commands/seurat-5.5.0/R/preprocessing.R:4865
LogNormalize.V3Matrix <- function(
  object,
  scale.factor = 1e4
){
    norm.data <- LogNorm(object, scale_factor = scale.factor)
    colnames(x = norm.data) <- colnames(x = object)
    rownames(x = norm.data) <- rownames(x = object)
    return(norm.data)
}

# commands/seurat-5.5.0/R/preprocessing5.R:311
NormalizeData.StdAssay <- function(
  object,
  scale.factor = 1e4
) {
    layer <- Layers(object = object, search = "counts")

    LayerData(
      object = object,
      layer = "data",
      features = Features(x = object, layer = "counts"),
      cells = Cells(x = object, layer = "counts")
    ) <- LogNormalize.V3Matrix(
      object = LayerData(object = object, layer = "counts", fast = NA),
      scale.factor = scale.factor
    )
    return(object)

}

# commands/seurat-5.5.0/R/preprocessing.R:5056
NormalizeData.Seurat <- function(
  object,
  scale.factor = 1e4
) {
    assay <- DefaultAssay(object = object)
    assay.data <- NormalizeData.StdAssay(
        object = object[[assay]],
        scale.factor = scale.factor
    )
    object[[assay]] <- assay.data
    return(object)
}

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Normalise
pbmc <- NormalizeData.Seurat(pbmc, scale.factor = 10000)


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
# TMSB4X         5.496190         4.825385         5.050981         5.251126
# MALAT1         3.930709         3.226986         3.970951         3.637138
# B2M            4.152408         4.549060         4.809943         4.743789
# RPL13A         4.668858         4.243997         4.449180         4.016405
# RPL10          4.308571         4.528085         4.296976         3.980704
#        GGCACGTGTGAGAA-1
# TMSB4X         5.345753
# MALAT1         3.936645
# B2M            4.348158
# RPL13A         4.591010
# RPL10          4.761774