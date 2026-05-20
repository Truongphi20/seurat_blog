library(Seurat)
library(ggplot2)
library(sctransform)

# commands/seurat-5.5.0/R/preprocessing.R:6105
RegressOutMatrix <- function(
  data.expr,
  latent.data = NULL,
  features.regress = NULL
) {
    features.regress <- intersect(x = features.regress, y = rownames(x = data.expr))

    # Create formula for regression
    vars.to.regress <- colnames(x = latent.data)
    fmla <- paste('GENE ~', paste(vars.to.regress, collapse = '+'))
    fmla <- as.formula(object = fmla)

    regression.mat <- cbind(latent.data, data.expr[1,])
    colnames(regression.mat) <- c(colnames(x = latent.data), "GENE")
    qr <- lm(fmla, data = regression.mat, qr = TRUE)$qr

    # Make results matrix
    data.resid <- matrix(
        nrow = nrow(x = data.expr),
        ncol = ncol(x = data.expr)
    )

    for (i in 1:length(x = features.regress)) {
        x <- features.regress[i]
        regression.mat <- cbind(latent.data, data.expr[x, ])
        colnames(x = regression.mat) <- c(vars.to.regress, 'GENE')
        regression.mat <- qr.resid(qr = qr, y = data.expr[x,])
        data.resid[i, ] <- regression.mat
    }

    dimnames(x = data.resid) <- dimnames(x = data.expr)
    return(data.resid)
}

# commands/seurat-5.5.0/R/preprocessing.R:5111
ScaleData.default <- function(
  object,
  latent.data = NULL,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000
){
    features <- rownames(x = object)
    features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
    object <- object[features, , drop = FALSE]
    object.names <- dimnames(x = object)
    min.cells.to.block <- min(min.cells.to.block, ncol(x = object))

    split.by <- TRUE
    split.cells <- split(x = colnames(x = object), f = split.by)
    latent.data <- latent.data[colnames(x = object), , drop = FALSE]
    rownames(x = latent.data) <- colnames(x = object)

    object <- lapply(
        X = names(x = split.cells),
        FUN = function(x) {
            return(RegressOutMatrix(
                data.expr = object[, split.cells[[x]], drop = FALSE],
                latent.data = latent.data[split.cells[[x]], , drop = FALSE],
                features.regress = features
            ))
        }
    )
    object <- do.call(what = 'cbind', args = object)
    dimnames(x = object) <- object.names

    object <- as.matrix(x = object)
    scale.function <- FastRowScale

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
    }

    dimnames(x = scaled.data) <- object.names
    scaled.data[is.na(x = scaled.data)] <- 0

    return(scaled.data)
}

# commands/seurat-5.5.0/R/preprocessing.R:3863
SCTransform.default <- function(
  object,
  cell.attr,
  do.correct.umi = TRUE,
  ncells = 5000,
  variable.features.n = 3000,
  latent.data = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = umi) / 30), sqrt(x = ncol(x = umi) / 30)),
  vst.flavor = 'v2',
  seed.use = 1448145
) {
    set.seed(seed = seed.use)
    vst.args <- list()
    object <- as.sparse(x = object)
    umi <- object

    vst.args[['vst.flavor']] <- vst.flavor
    vst.args[['umi']] <- umi
    vst.args[['cell_attr']] <- cell.attr
    vst.args[['verbosity']] <- as.numeric(x = TRUE) * 1
    vst.args[['return_cell_attr']] <- TRUE
    vst.args[['return_gene_attr']] <- TRUE
    vst.args[['return_corrected_umi']] <- do.correct.umi
    vst.args[['n_cells']] <- min(ncells, ncol(x = umi))

    vst.out <- do.call(what = 'vst', args = vst.args)

    feature.variance <- vst.out$gene_attr[,"residual_variance"]
    names(x = feature.variance) <- rownames(x = vst.out$gene_attr)

    feature.variance <- sort(x = feature.variance, decreasing = TRUE)
    top.features <- names(x = feature.variance)[1:min(variable.features.n, length(x = feature.variance))]
    vst.out$y <- vst.out$y[top.features, ]

    scale.data <- vst.out$y
    # clip the residuals
    scale.data[scale.data < clip.range[1]] <- clip.range[1]
    scale.data[scale.data > clip.range[2]] <- clip.range[2]
    # 2nd regression
    scale.data <- ScaleData.default(
        scale.data,
        latent.data = latent.data,
        do.scale = do.scale,
        do.center = do.center,
        scale.max = Inf,
        block.size = 750,
        min.cells.to.block = 3000
    )
    vst.out$y <- scale.data
    vst.out$variable_features <- top.features
    return(vst.out)
}

# commands/seurat-5.5.0/R/preprocessing.R:4136
SCTransform.Assay <- function(
    object,
    cell.attr,
    do.correct.umi = TRUE,
    ncells = 5000,
    variable.features.n = 3000,
    latent.data = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
    vst.flavor = 'v2',
    seed.use = 1448145
) {
    set.seed(seed = seed.use)
    umi <- GetAssayData(object = object, layer = 'counts')
    vst.out <- SCTransform.default(object = umi,
                         cell.attr = cell.attr,
                         do.correct.umi = do.correct.umi,
                         ncells = ncells,
                         variable.features.n = variable.features.n,
                         latent.data = latent.data,
                         do.scale = do.scale,
                         do.center = do.center,
                         clip.range = clip.range,
                         vst.flavor = vst.flavor,
                         seed.use = seed.use)
    
    sct.method = NULL

    # create output assay and put (corrected) umi counts in count slot
    assay.out <- CreateAssayObject(counts = vst.out$umi_corrected)
    vst.out$umi_corrected <- NULL

    # set the variable genes
    VariableFeatures(object = assay.out) <- vst.out$variable_features
    # put log1p transformed counts in data
    assay.out <- SetAssayData(
        object = assay.out,
        layer = 'data',
        new.data = log1p(x = GetAssayData(object = assay.out, layer = 'counts'))
    )
    scale.data <- vst.out$y
    assay.out <- SetAssayData(
        object = assay.out,
        layer = 'scale.data',
        new.data = scale.data
    )
    vst.out$y <- NULL
    # save clip.range into vst model
    vst.out$arguments$sct.clip.range <- clip.range
    vst.out$arguments$sct.method <- sct.method
    Misc(object = assay.out, slot = 'vst.out') <- vst.out
    assay.out <- as(object = assay.out, Class = "SCTAssay")
    return(assay.out)

}

# commands/seurat-5.5.0/R/preprocessing5.R:1115
SCTransform.StdAssay <- function(
  object,
  cell.attr = NULL,
  do.correct.umi = TRUE,
  ncells = 5000,
  variable.features.n = 3000,
  latent.data = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
  vst.flavor = 'v2',
  seed.use = 1448145
){
    # Extract TK and TK.
    layer_name = "counts"
    layer_counts <- LayerData(object, layer = layer_name)
    layer_object <- CreateAssayObject(layer_counts)


    # Apply SCTransform to each assay in `input_list`.
    .cell.attr <- cell.attr[Cells(layer_object), ]
    assay_out <- SCTransform.Assay(
        layer_object,
        cell.attr = .cell.attr,
        do.correct.umi = do.correct.umi,
        ncells = ncells,
        variable.features.n = variable.features.n,
        latent.data = latent.data,
        do.scale = do.scale,
        do.center = do.center,
        clip.range = clip.range,
        vst.flavor = vst.flavor,
        seed.use = seed.use
    )

    var_features_union <- VariableFeatures(assay_out)
    all_features_intersect <- rownames(assay_out)

    # Keep features that are variable in at least one output assay/layer but
    # present in all of them.
    scale_data_features <- intersect(all_features_intersect, var_features_union)

    # Extract residuals for the selected features and store them in
    # the outputs scaled.data slot.
    residuals <- FetchResiduals(
        object = assay_out, 
        umi.object = object,
        features = scale_data_features,
        verbose = FALSE
    )
    
    LayerData(assay_out, layer = "scale.data") <- residuals

    # Set the output's variable features.
    VariableFeatures(assay_out) <- VariableFeatures(
        assay_out, 
        use.var.features = FALSE,
        nfeatures = variable.features.n
    )

    return(assay_out)
}

# commands/seurat-5.5.0/R/preprocessing.R:4235
SCTransform.Seurat <- function(
    object,
    new.assay.name = 'SCT',
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    vars.to.regress = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object[[assay]]) / 30), sqrt(x = ncol(x = object[[assay]]) / 30)),
    vst.flavor = "v2",
    seed.use = 1448145
){
    set.seed(seed = seed.use)
    vars.to.regress.subset <- vars.to.regress[vars.to.regress %in% colnames(x = object[[]])]
    latent.data <- object[[vars.to.regress.subset]]

    assay = "RNA"
    cell.attr <- slot(object = object, name = 'meta.data')[colnames(object[[assay]]),]
    assay.data <- SCTransform.StdAssay(object = object[[assay]],
                                cell.attr = cell.attr,
                                do.correct.umi = do.correct.umi,
                                ncells = ncells,
                                variable.features.n = variable.features.n,
                                latent.data = latent.data,
                                do.scale = do.scale,
                                do.center = do.center,
                                clip.range = clip.range,
                                vst.flavor = vst.flavor,
                                seed.use = seed.use)

    object[[new.assay.name]] <- assay.data
    
    return(object)
}


pbmc_data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc_data)

# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
pbmc <- SCTransform.Seurat(pbmc, vars.to.regress = "percent.mt")


## Checking the result 
# Extract the full matrix
mat <- pbmc@assays$SCT@counts

# Find rows (genes) and columns (cells) with the highest expressions
top_genes <- order(rowSums(mat > 0), decreasing = TRUE)[1:5]
top_cells <- order(colSums(mat > 0), decreasing = TRUE)[1:5]

# Subset the matrix using these top indices
dense_chunk <- mat[top_genes, top_cells]

# Print it as a standard, readable matrix
print(as.matrix(dense_chunk))

# Expected output 
#        TAATGCCTCGTCTC-1 CACCGGGACTTCTA-1 ACGTCGCTCCTGAA-1 CAGTTTACACACGT-1
# RPL15                11                8                6               14
# RPS18                 8                9                2               19
# RPS12                 6                4                4               12
# TMSB4X               47               75               24               62
# RPL10                20               13               17               19
#        GGGAACGAAGCTCA-1
# RPL15                 5
# RPS18                 5
# RPS12                 6
# TMSB4X               68
# RPL10                14