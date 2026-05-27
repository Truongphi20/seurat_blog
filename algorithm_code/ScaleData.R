library(dplyr)
library(Seurat)
library(patchwork)


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
pbmc <- ScaleData(pbmc, features = all.genes)

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