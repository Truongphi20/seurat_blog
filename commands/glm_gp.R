library(glmGamPoi)

my_stack <- readRDS("/workspaces/seurat_blog/test_data/run_glm_gp.rds")

fit <- glmGamPoi::glm_gp(data = my_stack$umi,
                           design = as.formula(my_stack$new_formula),
                           col_data = my_stack$data,
                           offset = my_stack$log_umi,
                           size_factors = FALSE)

fit$theta <- 1 / fit$overdispersions

print(head(fit$theta))