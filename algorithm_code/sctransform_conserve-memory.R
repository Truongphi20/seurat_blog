library(Seurat)
library(SeuratObject)
library(ggplot2)
library(sctransform)

# commands/sctransform-0.4.3/R/utils.R:460
get_model_formula <- function(model_str) {
  as.formula(gsub('^y', '', model_str))
}

# commands/sctransform-0.4.3/R/utils.R:473
prepare_regressor_data <- function(vst_out, cell_attr) {
  regressor_data <- model.matrix(get_model_formula(vst_out$model_str), cell_attr)

  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(get_model_formula(vst_out$model_str_nonreg), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
  }

  return(regressor_data)
}

# commands/sctransform-0.4.3/R/utils.R:428
get_nz_median2 <- function(umi, genes = NULL){
  if (is.null(genes)) {
    # Compute median for the entire matrix
    return(median(umi@x))
  } else if (length(genes) == 1) {
    # If only one gene is being subsetted
    return(median(umi[genes, umi[genes,] != 0]))
  } else if (length(genes) > 1) {
    # If multiple genes are being subsetted
    return(median(umi[genes,]@x))
  } else {
    stop("genes does not contain a vector of gene names")
  }
}

# commands/sctransform-0.4.3/R/utils.R:473
prepare_regressor_data <- function(vst_out, cell_attr) {
  regressor_data <- model.matrix(get_model_formula(vst_out$model_str), cell_attr)

  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(get_model_formula(vst_out$model_str_nonreg), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
  }

  return(regressor_data)
}

# commands/sctransform-0.4.3/R/utils.R:491
setup_progress_bar <- function(n_items, bin_size) {
  bin_ind <- ceiling(x = 1:n_items / bin_size)
  max_bin <- max(bin_ind)
  pb <- NULL

  return(list(bin_ind = bin_ind, pb = pb, max_bin = max_bin))
}

# commands/sctransform-0.4.3/R/utils.R:166
pearson_residual <- function(y, mu, theta, min_var = -Inf) {
  model_var <- mu + mu^2 / theta
  model_var[model_var < min_var] <- min_var
  return((y - mu) / sqrt(model_var))
}

# commands/sctransform-0.4.3/R/utils.R:530
clip_matrix_values <- function(mat, clip_range) {
  mat[mat < clip_range[1]] <- clip_range[1]
  mat[mat > clip_range[2]] <- clip_range[2]
  return(mat)
}


# commands/sctransform-0.4.3/R/utils.R:215
get_residuals <- function(vst_out, umi, residual_type = 'pearson',
                          res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                          min_variance = vst_out$arguments$min_variance,
                          cell_attr = vst_out$cell_attr, bin_size = 256,
                          verbosity = vst_out$arguments$verbosity)
{
    # Maximum pearson residual for non-zero median UMI is 5
    min_var <- (get_nz_median2(umi) / 5)^2

    regressor_data <- prepare_regressor_data(vst_out, cell_attr)
    model_pars <- vst_out$model_pars_fit

    genes <- rownames(umi)[rownames(umi) %in% rownames(model_pars)]

    pb_setup <- setup_progress_bar(length(genes), bin_size)
    bin_ind <- pb_setup$bin_ind
    max_bin <- pb_setup$max_bin
    res <- matrix(NA_real_, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))

    for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))

    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    res[genes_bin, ] <- pearson_residual(y, mu, model_pars[genes_bin, 'theta'], min_var = min_var)
  }

  res <- clip_matrix_values(res, res_clip_range)
  return(res)
}

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
) {
    features <- rownames(x = object)
    features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
    object <- object[features, , drop = FALSE]
    object.names <- dimnames(x = object)
    min.cells.to.block <- min(min.cells.to.block, ncol(x = object))
    split.by <- TRUE
    split.cells <- split(x = colnames(x = object), f = split.by)

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
  reference.SCT.model = NULL,
  do.correct.umi = TRUE,
  ncells = 5000,
  residual.features = NULL,
  variable.features.n = 3000,
  variable.features.rv.th = 1.3,
  vars.to.regress = NULL,
  latent.data = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = umi) / 30), sqrt(x = ncol(x = umi) / 30)),
  vst.flavor = 'v2',
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = 1448145,
  verbose = TRUE,
  ...
){
    set.seed(seed = seed.use)
    vst.args <- list(...)
    object <- as.sparse(x = object)
    umi <- object

    vst.args[['vst.flavor']] <- vst.flavor
    vst.args[['umi']] <- umi
    vst.args[['cell_attr']] <- cell.attr
    vst.args[['verbosity']] <- as.numeric(x = verbose) * 1
    vst.args[['return_cell_attr']] <- TRUE
    vst.args[['return_gene_attr']] <- TRUE
    vst.args[['return_corrected_umi']] <- do.correct.umi
    vst.args[['n_cells']] <- min(ncells, ncol(x = umi))
    residual.type <- 'pearson'
    res.clip.range <- c(-sqrt(x = ncol(x = umi)), sqrt(x = ncol(x = umi)))

    return.only.var.genes <- TRUE
    vst.args[['residual_type']] <- 'none'
    vst.out <- do.call(what = 'vst', args = vst.args)
    feature.variance <- get_residual_var(
        vst_out = vst.out,
        umi = umi,
        residual_type = residual.type,
        res_clip_range = res.clip.range
    )
    vst.out$gene_attr$residual_variance <- NA_real_
    vst.out$gene_attr[names(x = feature.variance), 'residual_variance'] <- feature.variance

    feature.variance <- vst.out$gene_attr[,"residual_variance"]
    names(x = feature.variance) <- rownames(x = vst.out$gene_attr)

    feature.variance <- sort(x = feature.variance, decreasing = TRUE)

    top.features <- names(x = feature.variance)[1:min(variable.features.n, length(x = feature.variance))]

    vst.out$y <- get_residuals(
        vst_out = vst.out,
        umi = umi[top.features, ],
        residual_type = residual.type,
        res_clip_range = res.clip.range,
        verbosity = as.numeric(x = verbose)*2
    )
    vst.out$gene_attr$residual_mean <- NA_real_
    vst.out$gene_attr[top.features, "residual_mean"] = matrixStats:::rowMeans2(x =  vst.out$y)
    
    scale.data <- vst.out$y
    # clip the residuals
    scale.data[scale.data < clip.range[1]] <- clip.range[1]
    scale.data[scale.data > clip.range[2]] <- clip.range[2]
    # 2nd regression
    scale.data <- ScaleData.default(
        scale.data,
        features = NULL,
        vars.to.regress = vars.to.regress,
        latent.data = latent.data,
        model.use = 'linear',
        use.umi = FALSE,
        do.scale = do.scale,
        do.center = do.center,
        scale.max = Inf,
        block.size = 750,
        min.cells.to.block = 3000,
        verbose = verbose
    )
    vst.out$y <- scale.data
    vst.out$variable_features <- top.features
    vst.out$umi_corrected <- umi
    
    return(vst.out)
}

# commands/seurat-5.5.0/R/preprocessing.R:4136
SCTransform.Assay <- function(
    object,
    cell.attr,
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    vars.to.regress = NULL,
    latent.data = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
    vst.flavor = 'v2',
    conserve.memory = FALSE,
    return.only.var.genes = TRUE,
    seed.use = 1448145,
    verbose = TRUE,
    ...
){
    set.seed(seed = seed.use)
    do.correct.umi <- FALSE
    do.center <- FALSE

    umi <- GetAssayData(object = object, layer = 'counts')
    vst.out <- SCTransform.default(object = umi,
                         cell.attr = cell.attr,
                         reference.SCT.model = reference.SCT.model,
                         do.correct.umi = do.correct.umi,
                         ncells = ncells,
                         residual.features = residual.features,
                         variable.features.n = variable.features.n,
                         variable.features.rv.th = variable.features.rv.th,
                         vars.to.regress = vars.to.regress,
                         latent.data = latent.data,
                         do.scale = do.scale,
                         do.center = do.center,
                         clip.range = clip.range,
                         vst.flavor = vst.flavor,
                         conserve.memory = conserve.memory,
                         return.only.var.genes = return.only.var.genes,
                         seed.use = seed.use,
                         verbose = verbose,
                         ...)
        
  sct.method <- NULL

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
  layer = 'counts',
  cell.attr = NULL,
  reference.SCT.model = NULL,
  do.correct.umi = TRUE,
  ncells = 5000,
  residual.features = NULL,
  variable.features.n = 3000,
  variable.features.rv.th = 1.3,
  vars.to.regress = NULL,
  latent.data = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
  vst.flavor = 'v2',
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = 1448145,
  verbose = TRUE,
  ...
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
        reference.SCT.model = reference.SCT.model,
        do.correct.umi = do.correct.umi,
        ncells = ncells,
        residual.features = residual.features,
        variable.features.n = variable.features.n,
        variable.features.rv.th = variable.features.rv.th,
        vars.to.regress = vars.to.regress,
        latent.data = latent.data,
        do.scale = do.scale,
        do.center = do.center,
        clip.range = clip.range,
        vst.flavor = vst.flavor,
        conserve.memory = conserve.memory,
        return.only.var.genes = return.only.var.genes,
        seed.use = seed.use,
        verbose = verbose,
        ...
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
    assay = "RNA",
    new.assay.name = 'SCT',
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    vars.to.regress = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object[[assay]]) / 30), sqrt(x = ncol(x = object[[assay]]) / 30)),
    vst.flavor = "v2",
    conserve.memory = FALSE,
    return.only.var.genes = TRUE,
    seed.use = 1448145,
    verbose = TRUE,
    ...
){
    set.seed(seed = seed.use)
    latent.data <- NULL
    assay <- "RNA"

    cell.attr <- slot(object = object, name = 'meta.data')[colnames(object[[assay]]),]

    assay.data <- SCTransform.StdAssay(object = object[[assay]],
                        cell.attr = cell.attr,
                        reference.SCT.model = reference.SCT.model,
                        do.correct.umi = do.correct.umi,
                        ncells = ncells,
                        residual.features = residual.features,
                        variable.features.n = variable.features.n,
                        variable.features.rv.th = variable.features.rv.th,
                        vars.to.regress = vars.to.regress,
                        latent.data = latent.data,
                        do.scale = do.scale,
                        do.center = do.center,
                        clip.range = clip.range,
                        vst.flavor = vst.flavor,
                        conserve.memory = conserve.memory,
                        return.only.var.genes = return.only.var.genes,
                        seed.use = seed.use,
                        verbose = verbose,
                        ...)

  # Extract all SCT models stored in assay
  sct_models <- slot(object = assay.data, name = "SCTModel.list")
  
  # Update umi.assay field for every SCT model 
  slot(object = assay.data, name = "SCTModel.list") <- lapply(sct_models, function(model) {
    slot(model, name = "umi.assay") <- assay
    model
  })

  object[[new.assay.name]] <- assay.data
  return(object)
}

pbmc_data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc_data)

# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
pbmc <- SCTransform.Seurat(
  object = pbmc,
  method = 'glmGamPoi', 
  clip.range = c(-10000, 10),
  do.correct.umi = FALSE,
  conserve.memory = TRUE
)

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
#        CCAGTCTGCGGAGA-1 TTACTCGAACGTTG-1 AGAGGTCTACAGCT-1 GCGAAGGAGAGCTT-1
# TMSB4X              384              189              167              195
# MALAT1               79               37               56               38
# B2M                  99              143              131              117
# RPL13A              167              105               91               56
# RPL10               116              140               78               54
#        GGCACGTGTGAGAA-1
# TMSB4X              216
# MALAT1               52
# B2M                  79
# RPL13A              101
# RPL10               120