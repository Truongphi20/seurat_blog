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

# commands/sctransform-0.4.3/src/utils.cpp:88
row_gmean_dgcmatrix <- function(matrix, eps) {
  x <- matrix@x
  i <- matrix@i
  dim <- matrix@Dim
  rows <- dim[1]
  cols <- dim[2]
  
  # Initialize vectors
  ret <- numeric(rows) 
  nzero <- rep(cols, rows)
  
  log_eps <- log(eps)
  
  # Accumulate log values for non-zero elements
  log_x_plus_eps <- log(x + eps)
  
  agg_ret <- rowsum(log_x_plus_eps, i + 1)
  agg_nzero <- rowsum(rep(1, length(x)), i + 1)
  
  # Map the aggregated values back to their specific rows
  present_rows <- as.integer(rownames(agg_ret))
  ret[present_rows] <- agg_ret[, 1]
  nzero[present_rows] <- cols - agg_nzero[, 1]
  
  # Compute final geometric mean per row
  ret <- exp((ret + log_eps * nzero) / cols) - eps
  
  # Assign Row Names
  dn <- matrix@Dimnames
  if (!is.null(dn[[1]])) {
    names(ret) <- dn[[1]]
  }
  
  return(ret)
}

# commands/sctransform-0.4.3/R/utils.R:73
row_gmean <- function(x, eps = 1) {
    ret <- row_gmean_dgcmatrix(matrix = x, eps = eps)
    names(ret) <- rownames(x)
    return(ret)
}

# commands/seurat-5.5.0/src/stats.cpp:7
row_sum_dgcmatrix_r <- function(x, i, rows, cols) {
  rowsum_vec <- numeric(rows)
  
  # Accumulate values of non-zero elements into their respective rows
  agg_sum <- rowsum(x, i + 1)
  
  # Extract which rows actually had non-zero values
  present_rows <- as.integer(rownames(agg_sum))
  rowsum_vec[present_rows] <- agg_sum[, 1]
  
  return(rowsum_vec)
}

# commands/seurat-5.5.0/src/stats.cpp:18
row_mean_dgcmatrix_r <- function(x, i, rows, cols) {
  rowmean <- row_sum_dgcmatrix_r(x, i, rows, cols)
  
  return(rowmean / cols)
}

# commands/seurat-5.5.0/src/stats.cpp:28
row_var_dgcmatrix_r <- function(x, i, rows, cols) {
  # Get the row means using our pure R mean function
  rowmean <- row_mean_dgcmatrix_r(x, i, rows, cols)
  
  # Initialize output vectors
  rowvar <- numeric(rows)
  nzero <- rep(cols, rows)
  
  # Calculate squared differences for non-zero items
  sq_diff <- (x - rowmean[i + 1])^2
  
  agg_var <- rowsum(sq_diff, i + 1)
  agg_nzero <- rowsum(rep(1, length(x)), i + 1)
  
  present_rows <- as.integer(rownames(agg_var))
  rowvar[present_rows] <- agg_var[, 1]
  nzero[present_rows] <- cols - agg_nzero[, 1]
  
  # Add contribution of structural zeros and divide by (cols - 1)
  rowvar <- (rowvar + ((rowmean^2) * nzero)) / (cols - 1)
  
  return(rowvar)
}

# commands/sctransform-0.4.3/R/utils.R:94
row_var <- function(x) {
  if (inherits(x = x, what = 'matrix')) {
    ret <- matrixStats:::rowVars(x)
    names(ret) <- rownames(x)
    return(ret)
  }
  if (inherits(x = x, what = 'dgCMatrix')) {
    ret <- row_var_dgcmatrix_r(x = x@x, i = x@i, rows = nrow(x), cols = ncol(x))
    names(ret) <- rownames(x)
    return(ret)
  }
  stop('matrix x needs to be of class matrix or dgCMatrix')
}

# commands/sctransform-0.4.3/R/utils.R:3
make_cell_attr <- function(umi, cell_attr, latent_var, batch_var, latent_var_nonreg, verbosity) {

    new_attr <- list()
    new_attr$umi <- colSums(umi)
    new_attr$log_umi <- log10(new_attr$umi)

    new_attr <- do.call(cbind, new_attr)
    cell_attr <- cbind(cell_attr, new_attr[, setdiff(colnames(new_attr), colnames(cell_attr)), drop = FALSE])

    return(cell_attr)
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

# commands/sctransform-0.4.3/R/denoise.R:81
correct <- function(x, data = 'y', cell_attr = x$cell_attr, as_is = FALSE,
                    do_round = TRUE, do_pos = TRUE, scale_factor=NA, verbosity = 2)
{
    data <- x[[data]]
    latent_var <- x$arguments$latent_var

    cell_attr[, latent_var] <- apply(cell_attr[, latent_var, drop=FALSE], 2, function(x) rep(median(x), length(x)))

    regressor_data <- model.matrix(as.formula(gsub('^y', '',x$model_str)), cell_attr)
    genes <- rownames(data)
    bin_size <- x$arguments$bin_size

    bin_ind <- ceiling(x = 1:length(genes) / bin_size)
    max_bin <- max(bin_ind)
    corrected_data <- matrix(NA_real_, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))

    for (i in 1:max_bin) {
        genes_bin <- genes[bin_ind == i]
        pearson_residual <- data[genes_bin, ]
        coefs <- x$model_pars_fit[genes_bin, -1, drop=FALSE]
        theta <- x$model_pars_fit[genes_bin, 1]
        mu <- exp(tcrossprod(coefs, regressor_data))
        variance <- mu + mu^2 / theta
        corrected_data[genes_bin, ] <- mu + pearson_residual * sqrt(variance)
    }

    corrected_data <- round(corrected_data, 0)
    corrected_data[corrected_data < 0] <- 0

    return(corrected_data)
}

# commands/sctransform-0.4.3/R/vst.R:109
vst <- function(umi,
                cell_attr = NULL,
                latent_var = c('log_umi'),
                batch_var = NULL,
                latent_var_nonreg = NULL,
                n_genes = 2000,
                n_cells = NULL,
                method = 'poisson',
                do_regularize = TRUE,
                theta_regularization = 'od_factor',
                res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                bin_size = 500,
                min_cells = 5,
                residual_type = 'pearson',
                return_cell_attr = FALSE,
                return_gene_attr = TRUE,
                return_corrected_umi = FALSE,
                min_variance = -Inf,
                bw_adjust = 3,
                gmean_eps = 1,
                theta_estimation_fun = 'theta.ml',
                theta_given = NULL,
                exclude_poisson = FALSE,
                use_geometric_mean = TRUE,
                use_geometric_mean_offset = FALSE,
                fix_intercept = FALSE,
                fix_slope = FALSE,
                scale_factor = NA,
                vst.flavor = NULL,
                verbosity = 2)
{
    method <- "glmGamPoi_offset"
    exclude_poisson <- TRUE
    min_variance <- 'umi_median'

    arguments <- as.list(environment())
    arguments <- arguments[!names(arguments) %in% c("umi", "cell_attr")]

    cell_attr <- make_cell_attr(umi, cell_attr, latent_var, batch_var, latent_var_nonreg, verbosity)

    # we will generate output for all genes detected in at least min_cells cells
    # but for the first step of parameter estimation we might use only a subset of genes
    genes_cell_count <- rowSums(umi >= 0.01)
    genes <- rownames(umi)[genes_cell_count >= min_cells]
    umi <- umi[genes, ]

    genes_log_gmean <- log10(row_gmean(umi, eps = gmean_eps))

    cells_step1 <- colnames(umi)
    genes_step1 <- genes
    genes_log_gmean_step1 <- genes_log_gmean

    genes_amean <- NULL
    genes_var <- NULL
    # Exclude known poisson genes from the learning step
    genes_amean <- rowMeans(umi)
    genes_var <- row_var(umi)
    overdispersion_factor <- genes_var - genes_amean
    overdispersion_factor_step1 <- overdispersion_factor[genes_step1]
    is_overdispersed <- (overdispersion_factor_step1 > 0)

    genes_step1 <-  genes_step1[is_overdispersed]
    genes_log_gmean_step1 <-  genes_log_gmean[genes_step1]

    data_step1 <- cell_attr[cells_step1, , drop = FALSE]

    # density-sample genes to speed up the first step
    log_gmean_dens <- density(x = genes_log_gmean_step1, bw = 'nrd', adjust = 1)
    sampling_prob <- 1 / (approx(x = log_gmean_dens$x, y = log_gmean_dens$y, xout = genes_log_gmean_step1)$y + .Machine$double.eps)
    genes_step1 <- sample(x = genes_step1, size = n_genes, prob = sampling_prob)
    genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, ], eps = gmean_eps))

    model_str <- paste0('y ~ ', paste(latent_var, collapse = ' + '))

    bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
    max_bin <- max(bin_ind)

    model_pars <- sctransform:::get_model_pars(genes_step1, bin_size, umi, model_str, cells_step1,
                               method, data_step1, theta_given, theta_estimation_fun,
                               exclude_poisson, fix_intercept, fix_slope,
                               use_geometric_mean, use_geometric_mean_offset, verbosity)
    
    model_pars_fit <- sctransform:::reg_model_pars(model_pars, genes_log_gmean_step1, genes_log_gmean, cell_attr,
                                     batch_var, cells_step1, genes_step1, umi, bw_adjust, gmean_eps,
                                     theta_regularization, genes_amean, genes_var,
                                     exclude_poisson, fix_intercept, fix_slope,
                                     use_geometric_mean, use_geometric_mean_offset, verbosity)
    model_pars_outliers <- attr(model_pars_fit, 'outliers')

    # use all fitted values in NB model
    regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), cell_attr)

    model_str_nonreg <- ''
    model_pars_nonreg <- c()
    model_pars_final <- model_pars_fit
    regressor_data_final <- regressor_data

    # Maximum pearson residual for non-zero median UMI is 5
    min_var <- (get_nz_median2(umi) / 5)^2
    arguments$set_min_var <- min_var

    bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
    max_bin <- max(bin_ind)

    res <- matrix(NA_real_, length(genes), nrow(regressor_data_final), dimnames = list(genes, rownames(regressor_data_final)))
    for (i in 1:max_bin){
        genes_bin <- genes[bin_ind == i]
        mu <- exp(tcrossprod(model_pars_final[genes_bin, -1, drop=FALSE], regressor_data_final))
        y <- as.matrix(umi[genes_bin, , drop=FALSE])

        res[genes_bin, ] <- pearson_residual(y, mu, model_pars_final[genes_bin, 'theta'], min_var = min_var)
    }

    rv <- list(y = res,
             model_str = model_str,
             model_pars = model_pars,
             model_pars_outliers = model_pars_outliers,
             model_pars_fit = model_pars_fit,
             model_str_nonreg = model_str_nonreg,
             model_pars_nonreg = model_pars_nonreg,
             arguments = arguments,
             genes_log_gmean_step1 = genes_log_gmean_step1,
             cells_step1 = cells_step1,
             cell_attr = cell_attr)
    rm(res)

    rv$umi_corrected <- correct(rv, do_round = TRUE, do_pos = TRUE, scale_factor = scale_factor,
                                               verbosity = verbosity)

    rv$y <- clip_matrix_values(rv$y, res_clip_range)

    gene_attr <- data.frame(
      detection_rate = genes_cell_count[genes] / ncol(umi),
      gmean = 10 ^ genes_log_gmean,
      amean = rowMeans(umi),
      variance = row_var(umi)
    )
    
    if (ncol(rv$y) > 0) {
      gene_attr$residual_mean = rowMeans(rv$y)
      gene_attr$residual_variance = row_var(rv$y)
    }

    rv[['gene_attr']] <- gene_attr
    return(rv)
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