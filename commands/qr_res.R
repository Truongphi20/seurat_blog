# Load shared lib
stats_dll <- dyn.load("/workspaces/seurat_blog/commands/r-source-4.5.2/modules/lapack.so") # Build R first

# Load stack frame
my_stack <- readRDS("/workspaces/seurat_blog/test_data/inputs_qr-resid.rds")

# commands/r-source-4.5.2/src/library/base/R/qr.R:162
qr.resid_cloned <- function(qr, y)
{
    k <- as.integer(qr$rank)
    n <- as.integer(nrow(qr$qr))
    ny <- as.integer(NCOL(y))

    storage.mode(y) <- "double"
    .Fortran(stats_dll$dqrrsd, as.double(qr$qr), n, k, as.double(qr$qraux), y, ny, rsd = y)$rsd
}

# commands/seurat-5.5.0/R/preprocessing.R:6176
regression.mat = qr.resid_cloned(qr = my_stack$qr, y = my_stack$data.expr[my_stack$x,])

print(head(regression.mat))

## Expected output:
# AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1 
#        -1.962384        -2.289679        -1.493937         4.740084 
# AAACCGTGTATGCG-1 AAACGCACTGGTAC-1 
#        -1.312671        -1.605116