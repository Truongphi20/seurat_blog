# Load shared lib
stats_dll <- dyn.load("/workspaces/seurat_blog/commands/r-source-4.5.2/src/library/stats/src/build/stats.so")

## Read stack variables
mystack <- readRDS("/workspaces/seurat_blog/test_data/inputs_ksmooth.rds")

# commands/r-source-4.5.2/src/library/stats/R/ksmooth.R:19
cloned_ksmooth <- function(x, y, kernel = c("box", "normal"), bandwidth = 0.5,
           range.x = range(x), n.points = max(100L, length(x)), x.points)
{
    browser()
    ## box is [-0.5, 0.5]. normal is sd = 1.4826/4
    if(missing(y) || is.null(y))
        stop("numeric y must be supplied.\nFor density estimation use density()")

    kernel <- match.arg(kernel)
    krn <- switch(kernel, "box" = 1L, "normal" = 2L)

    x.points <-
    if(missing(x.points))
        seq.int(range.x[1L], range.x[2L], length.out = n.points)
    else { n.points <- length(x.points); sort(x.points) }

    ord <- order(x)
    .Call(stats_dll$ksmooth, x[ord], y[ord], x.points, krn, bandwidth)
}

# algorithm_code/sctransform_default.R:472
output <- cloned_ksmooth(
    x = mystack$genes_log_gmean_step1, 
    y = mystack$model_pars[, 'dispersion_par'],
    x.points = mystack$x_points, 
    bandwidth = mystack$bw, 
    kernel='normal'
)

print(head(output$y))

## Expected output:
# [1] 0.1276975 0.1276975 0.1276975 0.1276975 0.1276975 0.1276975