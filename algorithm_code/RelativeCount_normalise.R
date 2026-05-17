library(dplyr)
library(Seurat)
library(patchwork)

NormalizeData.Seurat <- function(
  object,
  assay = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  margin = 1,
  verbose = TRUE,
  ...
){
    assay <- DefaultAssay(object = object)
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

# Normalise
pbmc <- NormalizeData.Seurat(pbmc, normalization.method = "RC", scale.factor = 10000)

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
# TMSB4X        242.76141        123.63446        155.17562        189.78102
# MALAT1         49.94310         24.20357         52.03494         36.98297
# B2M            62.58693         93.54353        121.72459        113.86861
# RPL13A        105.57593         68.68581         84.55677         54.50122
# RPL10          73.33418         91.58108         72.47723         52.55474
#        GGCACGTGTGAGAA-1
# TMSB4X        208.71582
# MALAT1         50.24640
# B2M            76.33588
# RPL13A         97.59397
# RPL10         115.95323