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

# commands/sctransform-0.4.3/R/utils.R:3
make_cell_attr <- function(umi, cell_attr, latent_var, batch_var, latent_var_nonreg, verbosity) {

    new_attr <- list()
    new_attr$umi <- colSums(umi)
    new_attr$log_umi <- log10(new_attr$umi)

    new_attr <- do.call(cbind, new_attr)
    cell_attr <- cbind(cell_attr, new_attr[, setdiff(colnames(new_attr), colnames(cell_attr)), drop = FALSE])

    return(cell_attr)
}

# commands/sctransform-0.4.3/R/utils.R:73
row_gmean <- function(x, eps = 1) {
    ret <- sctransform:::row_gmean_dgcmatrix(matrix = x, eps = eps)
    names(ret) <- rownames(x)
    return(ret)
}

# commands/sctransform-0.4.3/R/utils.R:94
row_var <- function(x) {
    ret <- sctransform:::row_var_dgcmatrix(x = x@x, i = x@i, rows = nrow(x), cols = ncol(x))
    names(ret) <- rownames(x)
    return(ret)
}

# commands/sctransform-0.4.3/R/vst.R:467
get_model_pars <- function(genes_step1, bin_size, umi, model_str, cells_step1,
                           method, data_step1, theta_given, theta_estimation_fun,
                           exclude_poisson = FALSE, fix_intercept = FALSE,
                           fix_slope = FALSE, use_geometric_mean = TRUE,
                           use_geometric_mean_offset = FALSE, verbosity = 0)
{
  bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  model_pars <- list()
  for (i in 1:max_bin) {
    genes_bin_regress <- genes_step1[bin_ind == i]
    umi_bin <- as.matrix(umi[genes_bin_regress, cells_step1, drop=FALSE])

    n_workers <- 1
    genes_per_worker <- nrow(umi_bin) / n_workers + .Machine$double.eps
    index_vec <- 1:nrow(umi_bin)
    index_lst <- split(index_vec, ceiling(index_vec/genes_per_worker))

    par_lst <- list()
    for (indices in index_lst){

      umi_bin_worker <- umi_bin[indices, , drop = FALSE]
      res <- sctransform:::fit_glmGamPoi_offset(umi = umi_bin_worker, model_str = model_str,
                                        data = data_step1, allow_inf_theta = exclude_poisson)
      par_lst[[length(par_lst) + 1]] <- res
    }

    model_pars[[i]] <- do.call(rbind, par_lst)
  }
  model_pars <- do.call(rbind, model_pars)

  rownames(model_pars) <- genes_step1
  colnames(model_pars)[1] <- 'theta'

  genes_amean <- rowMeans(umi)
  genes_var <- row_var(umi)

  genes_amean_step1 <- genes_amean[genes_step1]
  genes_var_step1 <- genes_var[genes_step1]

  predicted_theta <- genes_amean_step1^2/(genes_var_step1-genes_amean_step1)
  actual_theta <- model_pars[genes_step1, "theta"]
  diff_theta <- predicted_theta/actual_theta
  model_pars <- cbind(model_pars, diff_theta)

  # if the naive and estimated MLE are 1000x apart, set theta estimate to Inf
  diff_theta_index <- rownames(model_pars[model_pars[genes_step1, "diff_theta"]< 1e-3,])

  # Replace theta by infinity
  model_pars[diff_theta_index, 1] <- Inf
  # drop diff_theta column
  model_pars <- model_pars[, -dim(model_pars)[2]]

  return(model_pars)
}

# commands/sctransform-0.4.3/R/vst.R:710
reg_model_pars <- function(model_pars, genes_log_gmean_step1, genes_log_gmean, cell_attr,
                           batch_var, cells_step1, genes_step1, umi, bw_adjust, gmean_eps,
                           theta_regularization,
                           genes_amean = NULL, genes_var = NULL, exclude_poisson = FALSE,
                           fix_intercept = FALSE, fix_slope = FALSE, use_geometric_mean = TRUE,
                           use_geometric_mean_offset = FALSE, verbosity = 0) 
{
  genes <- names(genes_log_gmean)

  overdispersion_factor <- genes_var - genes_amean
  overdispersion_factor_step1 <- overdispersion_factor[genes_step1]

  all_poisson_genes <- genes[overdispersion_factor<=0]

  # also set genes with mean < 1e-3 as poisson
  low_mean_genes <- genes[genes_amean<1e-3]
  all_poisson_genes <- union(all_poisson_genes, low_mean_genes)


  poisson_genes_step1 <- genes_step1[overdispersion_factor_step1<=0]

  poisson_genes2 <- rownames(model_pars[!is.finite(model_pars[, 'theta']),])
  poisson_genes3 <- intersect(low_mean_genes, genes_step1)
  poisson_genes_step1 <- union(union(poisson_genes_step1, poisson_genes2),poisson_genes3)

  # Call offset model with theta=inf
  # only the slope and intercept are used downstream
  mean_cell_sum <- mean(x = colSums(umi))
  vst_out_offset <- cbind(rep(Inf, length(all_poisson_genes)),
                              log(genes_amean[all_poisson_genes]) - log(mean_cell_sum),
                              rep(log(10), length(all_poisson_genes) ))
  dimnames(vst_out_offset) <- list(all_poisson_genes, c('theta', '(Intercept)', 'log_umi'))

  dispersion_par <- rep(0, dim(vst_out_offset)[1])
  vst_out_offset <- cbind(vst_out_offset, dispersion_par)

  # we don't regularize theta directly
  # prior to v0.3 we regularized log10(theta)
  # now we transform to overdispersion factor
  # variance of NB is mu * (1 + mu / theta)
  # (1 + mu / theta) is what we call overdispersion factor here
  dispersion_par <- switch(theta_regularization,
                           'log_theta' = log10(model_pars[, 'theta']),
                           'od_factor' = log10(1 + 10^genes_log_gmean_step1 / model_pars[, 'theta']),
                           stop('theta_regularization ', theta_regularization, ' unknown - only log_theta and od_factor supported at the moment')
  )

  model_pars_all <- model_pars

  model_pars <- model_pars[, colnames(model_pars) != 'theta']
  model_pars <- cbind(dispersion_par, model_pars)

  # look for outliers in the parameters
  # outliers are those that do not fit the overall relationship with the mean at all
  outliers <- apply(model_pars, 2, function(y) sctransform:::is_outlier(y, genes_log_gmean_step1))
  outliers <- apply(outliers, 1, any)

  # also call theta=inf as outliers
  is_theta_inf <- !is.finite(model_pars_all[, "theta"])
  outliers <- outliers | is_theta_inf

  model_pars <- model_pars[!outliers, ]
  genes_step1 <- rownames(model_pars)
  genes_log_gmean_step1 <- genes_log_gmean_step1[!outliers]

  overdispersed_genes <- setdiff(rownames(model_pars), all_poisson_genes)
  model_pars <- model_pars[overdispersed_genes, ]
  genes_step1 <- rownames(model_pars)
  genes_log_gmean_step1 <- genes_log_gmean_step1[overdispersed_genes]

  # select bandwidth to be used for smoothing
  bw <- bw.SJ(genes_log_gmean_step1) * bw_adjust

  # for parameter predictions
  x_points <- pmax(genes_log_gmean, min(genes_log_gmean_step1))
  x_points <- pmin(x_points, max(genes_log_gmean_step1))

  # take results from step 1 and fit/predict parameters to all genes
  o <- order(x_points)
  model_pars_fit <- matrix(NA_real_, length(genes), ncol(model_pars),
                           dimnames = list(genes, colnames(model_pars)))

  # fit / regularize dispersion parameter
  model_pars_fit[o, 'dispersion_par'] <- ksmooth(x = genes_log_gmean_step1, y = model_pars[, 'dispersion_par'],
                                                 x.points = x_points, bandwidth = bw, kernel='normal')$y

  # global fit / regularization for all coefficients
  for (i in 2:ncol(model_pars)) {
    model_pars_fit[o, i] <- ksmooth(x = genes_log_gmean_step1, y = model_pars[, i],
                                    x.points = x_points, bandwidth = bw, kernel='normal')$y
  }

  dispersion_par <- rep(0, length(all_poisson_genes))
  model_pars_fit[all_poisson_genes, "dispersion_par"] <- dispersion_par

  theta <- 10^genes_log_gmean / (10^model_pars_fit[, 'dispersion_par'] - 1)

  model_pars_fit <- model_pars_fit[, colnames(model_pars_fit) != 'dispersion_par']
  model_pars_fit <- cbind(theta, model_pars_fit)

  for (col in intersect(colnames(x = model_pars_fit), colnames(x = vst_out_offset)) ){
    stopifnot(col %in% colnames(vst_out_offset))
    model_pars_fit[all_poisson_genes, col] <- vst_out_offset[all_poisson_genes, col]
  }

  attr(model_pars_fit, 'outliers') <- outliers
  return(model_pars_fit)
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

    model_pars <- get_model_pars(genes_step1, bin_size, umi, model_str, cells_step1,
                               method, data_step1, theta_given, theta_estimation_fun,
                               exclude_poisson, fix_intercept, fix_slope,
                               use_geometric_mean, use_geometric_mean_offset, verbosity)
    
    model_pars_fit <- reg_model_pars(model_pars, genes_log_gmean_step1, genes_log_gmean, cell_attr,
                                     batch_var, cells_step1, genes_step1, umi, bw_adjust, gmean_eps,
                                     theta_regularization, genes_amean, genes_var,
                                     exclude_poisson, fix_intercept, fix_slope,
                                     use_geometric_mean, use_geometric_mean_offset, verbosity)
    model_pars_outliers <- attr(model_pars_fit, 'outliers')

    model_str_nonreg <- ''
    model_pars_nonreg <- c()

    res <- matrix(data = NA, nrow = 0, ncol = 0)
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

    rv$y <- clip_matrix_values(rv$y, res_clip_range)
    gene_attr <- data.frame(
      detection_rate = genes_cell_count[genes] / ncol(umi),
      gmean = 10 ^ genes_log_gmean,
      amean = rowMeans(umi),
      variance = row_var(umi)
    )
    rv[['gene_attr']] <- gene_attr

    return(rv)
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