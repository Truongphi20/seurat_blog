library(glmGamPoi)
library(beachmat)

## Load unoptimized shared lib
so_path <- "/usr/local/lib/R/site-library/glmGamPoi/libs/glmGamPoi.so"
new_so_path <- "/workspaces/seurat_blog/commands/glmGamPoi_1.24.0/src/build/glmGamPoi.so"

if (as.character(getLoadedDLLs()[["glmGamPoi"]])[2] == so_path ) {
    dyn.unload(so_path)
}

# Load your optimized version
dyn.load(new_so_path)

## Read stack variables
mystack <- readRDS("/workspaces/seurat_blog/test_data/inputs_conventional_score_function_fast.rds")

## Run conventional_score_function_fast
conventional_score_function_fast <- function(y, mu, log_theta, model_matrix, do_cr_adj, unique_counts = as.numeric( c()), count_frequencies = as.numeric( c())) {
    .Call("_glmGamPoi_conventional_score_function_fast", y, mu, log_theta, model_matrix, do_cr_adj, unique_counts, count_frequencies)
}

# commands/glmGamPoi_1.24.0/R/overdispersion.R:190
far_left_value <- conventional_score_function_fast(
    mystack$y, 
    mu = mystack$mean_vector, 
    log_theta = log(1e-8),
    model_matrix = mystack$model_matrix, 
    do_cr_adj = mystack$do_cox_reid_adjustment, 
    mystack$tab[[1]], 
    mystack$tab[[2]]
)