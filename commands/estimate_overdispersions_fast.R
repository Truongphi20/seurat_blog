library(glmGamPoi)
library(beachmat)

## Load unoptimized shared lib
so_path <- "/usr/local/lib/R/site-library/glmGamPoi/libs/glmGamPoi.so"
if (is.loaded(so_path)) {
    dyn.unload(so_path)
}
dyn.load("/workspaces/seurat_blog/commands/glmGamPoi_1.24.0/src/build/glmGamPoi.so")

## Read stack variables
mystack <- readRDS("/workspaces/seurat_blog/test_data/inputs_estimate_overdispersions_fast.rds")

## Run estimation of overdispersions
## Copy from commands/glmGamPoi_1.24.0/R/overdispersion.R:114
est <- glmGamPoi:::estimate_overdispersions_fast(
    initializeCpp(mystack$y), 
    initializeCpp(mystack$mean), 
    mystack$model_matrix, 
    mystack$do_cox_reid_adjustment, 
    mystack$n_subsamples, 
    mystack$max_iter
)

print(head(est$estimate))