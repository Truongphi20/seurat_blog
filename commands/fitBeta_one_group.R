library(glmGamPoi)
library(beachmat)

## Load unoptimized shared lib
so_path <- "/usr/local/lib/R/site-library/glmGamPoi/libs/glmGamPoi.so"
if (is.loaded(so_path)) {
    dyn.unload(so_path)
}
dyn.load("/workspaces/seurat_blog/commands/glmGamPoi_1.24.0/src/build/glmGamPoi.so")

## Read stack variables
mystack <- readRDS("/workspaces/seurat_blog/test_data/inputs_fitBeta_one_group.rds")

## Run fitBeta
fitBeta_one_group <- function(Y, offset_matrix, thetas, beta_start_values, tolerance, maxIter) {
    .Call("_glmGamPoi_fitBeta_one_group", Y, offset_matrix, thetas, beta_start_values, tolerance, maxIter)
}

fitBeta_one_group(initializeCpp(mystack$Y),
                initializeCpp(mystack$offset), thetas = mystack$dispersions,
                beta_start_values = mystack$beta_group_init,
                tolerance = 1e-8, maxIter = 100)