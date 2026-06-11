library(Rcpp)

# Load shared lib
dyn.load("/workspaces/seurat_blog/commands/glmGamPoi_1.24.0/src/build/glmGamPoi.so")

## Read stack variables
mystack <- readRDS("/workspaces/seurat_blog/test_data/inputs_nlminb.rds")

## Declare function calls
conventional_loglikelihood_fast <- function(y, mu, log_theta, model_matrix, do_cr_adj, unique_counts = as.numeric( c()), count_frequencies = as.numeric( c())) {
    .Call("_glmGamPoi_conventional_loglikelihood_fast", y, mu, log_theta, model_matrix, do_cr_adj, unique_counts, count_frequencies)
}

conventional_score_function_fast <- function(y, mu, log_theta, model_matrix, do_cr_adj, unique_counts = as.numeric( c()), count_frequencies = as.numeric( c())) {
    .Call("_glmGamPoi_conventional_score_function_fast", y, mu, log_theta, model_matrix, do_cr_adj, unique_counts, count_frequencies)
}

conventional_deriv_score_function_fast <- function(y, mu, log_theta, model_matrix, do_cr_adj, unique_counts = as.numeric( c()), count_frequencies = as.numeric( c())) {
    .Call("_glmGamPoi_conventional_deriv_score_function_fast", y, mu, log_theta, model_matrix, do_cr_adj, unique_counts, count_frequencies)
}

## Run nlminb
## Copy from commands/glmGamPoi_1.24.0/R/overdispersion.R:205
res <- nlminb(start = log(mystack$start_value),
         objective = function(log_theta){
           - conventional_loglikelihood_fast(mystack$y, mu = mystack$mean_vector, log_theta = log_theta,
                             model_matrix = mystack$model_matrix, do_cr_adj = mystack$do_cox_reid_adjustment, mystack$tab[[1]], mystack$tab[[2]])
         }, gradient = function(log_theta){
           - conventional_score_function_fast(mystack$y, mu = mystack$mean_vector, log_theta = log_theta,
                             model_matrix = mystack$model_matrix, do_cr_adj = mystack$do_cox_reid_adjustment, mystack$tab[[1]], mystack$tab[[2]])
         }, hessian = function(log_theta){
           res <- conventional_deriv_score_function_fast(mystack$y, mu = mystack$mean_vector, log_theta = log_theta,
                             model_matrix = mystack$model_matrix, do_cr_adj = mystack$do_cox_reid_adjustment, mystack$tab[[1]], mystack$tab[[2]])
           matrix(- res, nrow = 1, ncol = 1)
         }, lower = log(1e-16), upper = log(1e16),
         control = list(iter.max = mystack$max_iter))

## Print result
print(res$par)