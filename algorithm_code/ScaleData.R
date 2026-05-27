library(dplyr)
library(Seurat)
library(patchwork)

# commands/seurat-5.5.0/R/preprocessing.R:5111
ScaleData.default <- function(
  object,
  features = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
){
    features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
    object <- object[features, , drop = FALSE]
    object.names <- dimnames(x = object)
    min.cells.to.block <- min(min.cells.to.block, ncol(x = object))

    split.by <- TRUE
    split.cells <- split(x = colnames(x = object), f = split.by)

    scale.function <- Seurat:::FastSparseRowScale

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
        arg.list <- list(
            mat = object[features[my.inds], split.cells[["TRUE"]], drop = FALSE],
            scale = do.scale,
            center = do.center,
            scale_max = scale.max,
            display_progress = FALSE
        )
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[my.inds], split.cells[["TRUE"]]])
        scaled.data[features[my.inds], split.cells[["TRUE"]]] <- data.scale
        rm(data.scale)
    }

    dimnames(x = scaled.data) <- object.names
    scaled.data[is.na(x = scaled.data)] <- 0

    return(scaled.data)
}

# commands/seurat-5.5.0/R/preprocessing5.R:359
ScaleData.StdAssay <- function(
  object,
  features = NULL,
  layer = 'data',
  vars.to.regress = NULL,
  latent.data = NULL,
  by.layer = FALSE,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale= TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  save = 'scale.data',
  verbose = TRUE,
  ...
){
    use.umi <- FALSE
    layer <- 'data'

    ldata <- LayerData(object = object, layer = layer, features = features)

    ldata <- ScaleData.default(
      object = ldata,
      features = features,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      split.by = split.by,
      model.use = model.use,
      use.umi = use.umi,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = scale.max,
      block.size = block.size,
      min.cells.to.block = min.cells.to.block,
      verbose = verbose,
      ...
    )

    LayerData(object = object, layer = save, features = rownames(ldata)) <- ldata

    return(object)
}

# commands/seurat-5.5.0/R/preprocessing.R:5520
ScaleData.Seurat <- function(
  object,
  features = NULL,
  assay = NULL,
  vars.to.regress = NULL,
  split.by = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
    assay <- "RNA"
    latent.data <- NULL

    assay.data <- ScaleData.StdAssay(
        # object = assay.data,
        object = object[[assay]],
        features = features,
        vars.to.regress = vars.to.regress,
        latent.data = latent.data,
        split.by = split.by,
        model.use = model.use,
        use.umi = use.umi,
        do.scale = do.scale,
        do.center = do.center,
        scale.max = scale.max,
        block.size = block.size,
        min.cells.to.block = min.cells.to.block,
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