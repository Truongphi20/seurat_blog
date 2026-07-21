# Load stack frame
my_stack <- readRDS("/workspaces/seurat_blog/test_data/inputs_lm_qr.rds")

# commands/r-source-4.5.2/src/library/stats/all.R:14975
lm_cloned <- function (formula, data, subset, weights, na.action,
        method = "qr", model = TRUE, x = FALSE, y = FALSE,
        qr = TRUE, singular.ok = TRUE, contrasts = NULL,
        offset, ...)
{
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    ## need stats:: for non-standard evaluation
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (method == "model.frame")
    return(mf)
    else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'", method),
                domain = NA)
    mt <- attr(mf, "terms") # allow model.frame to update it
    y <- model.response(mf, "numeric")
    ## avoid any problems with 1D or nx1 arrays by as.vector.
    w <- as.vector(model.weights(mf))
    if(!is.null(w) && !is.numeric(w))
        stop("'weights' must be a numeric vector")
    offset <- model.offset(mf)
    mlm <- is.matrix(y)
    ny <- if(mlm) nrow(y) else length(y)
    if(!is.null(offset)) {
        if(!mlm) offset <- as.vector(offset)
        if(NROW(offset) != ny)
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                          NROW(offset), ny), domain = NA)
    }

    if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if(mlm) matrix(NA_real_, 0, ncol(y))
                 else numeric(),
          residuals = y,
          fitted.values = 0 * y, weights = w, rank = 0L,
          df.residual = if(!is.null(w)) sum(w != 0) else ny)
        if(!is.null(offset)) {
            z$fitted.values <- offset
            z$residuals <- y - offset
        }
    }
    else {
        x <- model.matrix(mt, mf, contrasts)
        z <- if(is.null(w)) lm.fit(x, y, offset = offset,
                                   singular.ok=singular.ok, ...)
            else lm.wfit(x, y, w, offset = offset, singular.ok=singular.ok, ...)
    }
    class(z) <- c(if(mlm) "mlm", "lm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if (model)
        z$model <- mf
    if (ret.x)
        z$x <- x
    if (ret.y)
        z$y <- y
    if (!qr) z$qr <- NULL
        z
}



# Note: regression.mat was cut down to first 1000 entries 
# commands/seurat-5.5.0/R/preprocessing.R:6159
qr <- lm_cloned(my_stack$fmla, data = my_stack$regression.mat, qr = TRUE)$qr

print(head(qr$qr))

# Expected output
#                   (Intercept)   percent.mt
# AAACATACAACCAC-1 -31.62277660 -69.86318083
# AAACATTGAGCTAC-1   0.03162278 -30.93523258
# AAACATTGATCAGC-1   0.03162278  -0.04347159
# AAACCGTGCTTCCG-1   0.03162278  -0.01586833
# AAACCGTGTATGCG-1   0.03162278  -0.03267256
# AAACGCACTGGTAC-1   0.03162278  -0.01843817
