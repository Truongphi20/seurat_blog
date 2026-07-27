library("irlba")

# Load stack frame
my_stack <- readRDS("/workspaces/seurat_blog/test_data/inputs_irlba.rds")

## Note that values in object (genes x cells) are scaled values

## Run pca
pca.results <- irlba(A = t(x = my_stack$object), nv = my_stack$npcs)

## View output
print(dim(pca.results$v))
print(pca.results$v[1:5,1:5])

### Expected outputs (Note: the signs can be changed)
# 2000   50

#              [,1]        [,2]        [,3]        [,4]        [,5]
# [1,] -0.011726988  0.01523301 -0.14757216  0.11112354 0.001589976
# [2,] -0.118696692  0.02094446 -0.02730169 -0.06085506 0.057716238
# [3,]  0.008678533  0.05451734  0.05480080  0.06629901 0.008642370
# [4,] -0.119595112  0.01637529 -0.01491376 -0.04505086 0.030425229
# [5,]  0.015225200 -0.13030508  0.03287074  0.05735748 0.125760240 