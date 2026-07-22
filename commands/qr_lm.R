# Load shared lib
stats_dll <- dyn.load("/workspaces/seurat_blog/commands/r-source-4.5.2/library/stats/libs/stats.so")

# Load stack frame
my_stack <- readRDS("/workspaces/seurat_blog/test_data/inputs_lm_qr.rds")

# commands/r-source-4.5.2/src/library/stats/all.R:15048
lm.fit_cloned <- function(x, y, tol = 1e-07)
{
    n <- nrow(x)
    p <- ncol(x)
    z <- .Call(stats_dll$Cdqrls, x, y, tol, FALSE)
    coef <- z$coefficients
    pivot <- z$pivot

    ## careful here: the rank might be 0
    r1 <- seq_len(z$rank)
    dn <- colnames(x) %||% paste0("x", 1L:p)
    nmeffects <- c(dn[pivot[r1]], rep.int("", n - z$rank))
    r2 <- if(z$rank < p) (z$rank+1L):p else integer()

    coef[r2] <- NA
    ## avoid copy
    if(z$pivoted) coef[pivot] <- coef
    names(coef) <- dn
    names(z$effects) <- nmeffects
    
    z$coefficients <- coef
    r1 <- y - z$residuals

    ## avoid unnecessary copy
    if(z$pivoted) colnames(z$qr) <- colnames(x)[z$pivot]
    qr <- z[c("qr", "qraux", "pivot", "tol", "rank")]

    result <- c(z[c("coefficients", "residuals", "effects", "rank")],
      list(fitted.values = r1, assign = attr(x, "assign"),
       qr = structure(qr, class="qr"),
       df.residual = n - z$rank))

    return(result)
}

# commands/r-source-4.5.2/src/library/stats/all.R:14975
lm_cloned <- function(formula, data, model = TRUE)
{
    mf <- match.call()
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    ## need stats:: for non-standard evaluation
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms") # allow model.frame to update it
    y <- model.response(mf, "numeric")
    ## avoid any problems with 1D or nx1 arrays by as.vector.
    offset <- model.offset(mf)
    mlm <- is.matrix(y)

    x <- model.matrix(mt, mf, NULL)
    z <- lm.fit_cloned(x, y)
    
    class(z) <- c(if(mlm) "mlm", "lm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$terms <- mt
    z$model <- mf
    
    return(z)
}



# Note: regression.mat was cut down to first 1000 entries 
# commands/seurat-5.5.0/R/preprocessing.R:6159
qr <- lm_cloned(my_stack$fmla, data = my_stack$regression.mat)$qr

print(head(qr$qr))

# Expected output
#                   (Intercept)   percent.mt
# AAACATACAACCAC-1 -31.62277660 -69.86318083
# AAACATTGAGCTAC-1   0.03162278 -30.93523258
# AAACATTGATCAGC-1   0.03162278  -0.04347159
# AAACCGTGCTTCCG-1   0.03162278  -0.01586833
# AAACCGTGTATGCG-1   0.03162278  -0.03267256
# AAACGCACTGGTAC-1   0.03162278  -0.01843817
