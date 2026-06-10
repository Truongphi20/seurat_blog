library(beachmat)

# Load shared lib
dyn.load("/workspaces/seurat_blog/commands/glmGamPoi_1.24.0/src/build/glmGamPoi.so")

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

print(far_left_value)