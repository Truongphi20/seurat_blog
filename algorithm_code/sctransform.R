library(Seurat)
library(ggplot2)
library(sctransform)

# commands/seurat-5.5.0/R/preprocessing5.R:1115
SCTransform.StdAssay <- function(
  object,
  layer = 'counts',
  cell.attr = NULL,
  reference.SCT.model = NULL,
  do.correct.umi = TRUE,
  ncells = 5000,
  residual.features = NULL,
  variable.features.n = 3000,
  variable.features.rv.th = 1.3,
  vars.to.regress = NULL,
  latent.data = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-sqrt(x = ncol(x = object) / 30), sqrt(x = ncol(x = object) / 30)),
  vst.flavor = 'v2',
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = 1448145,
  verbose = TRUE,
  ...
){
    # Extract TK and TK.
    layer_names = "counts"
    input_list <- lapply(
        layer_names,
        function(layer_name) {
            layer_counts <- LayerData(object, layer = layer_name)
            layer_object <- CreateAssayObject(layer_counts)
            return(layer_object)
        }
    )

    # Apply SCTransform to each assay in `input_list`.
    output_list <- lapply(
        input_list,
        function(input) {
            .cell.attr <- cell.attr[Cells(input), ]
            result <- SCTransform(
                input,
                cell.attr = .cell.attr,
                reference.SCT.model = reference.SCT.model,
                do.correct.umi = do.correct.umi,
                ncells = ncells,
                residual.features = residual.features,
                variable.features.n = variable.features.n,
                variable.features.rv.th = variable.features.rv.th,
                vars.to.regress = vars.to.regress,
                latent.data = latent.data,
                do.scale = do.scale,
                do.center = do.center,
                clip.range = clip.range,
                vst.flavor = vst.flavor,
                conserve.memory = conserve.memory,
                return.only.var.genes = return.only.var.genes,
                seed.use = seed.use,
                verbose = verbose,
                ...
            )
        }
    )

    assay_out <- output_list[[1]]
    # Take the union of variable features across all output assays/layers.
    var_features_union <- Reduce(
      union,
      lapply(
        output_list,
        function(output) {
          return(VariableFeatures(output))
        }
      )
    )
    # Take the intersection of all features across all output assays/layers.
    all_features_intersect <- Reduce(
      intersect,
      lapply(
        output_list,
        function(output) {
          return(rownames(output))
        }
      )
    )

    # Keep features that are variable in at least one output assay/layer but
    # present in all of them.
    scale_data_features <- intersect(all_features_intersect, var_features_union)

    # Extract residuals for the selected features and store them in
    # the outputs scaled.data slot.
    residuals <- suppressWarnings(
        FetchResiduals(
        object = assay_out, 
        umi.object = object,
        features = scale_data_features,
        verbose = FALSE
        )
    )
    LayerData(assay_out, layer = "scale.data") <- residuals

    # Set the output's variable features.
    VariableFeatures(assay_out) <- VariableFeatures(
        assay_out, 
        use.var.features = FALSE,
        nfeatures = variable.features.n
    )

    return(assay_out)
}

# commands/seurat-5.5.0/R/preprocessing.R:4235
SCTransform.Seurat <- function(
    object,
    assay = "RNA",
    new.assay.name = 'SCT',
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = 5000,
    residual.features = NULL,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    vars.to.regress = NULL,
    do.scale = FALSE,
    do.center = TRUE,
    clip.range = c(-sqrt(x = ncol(x = object[[assay]]) / 30), sqrt(x = ncol(x = object[[assay]]) / 30)),
    vst.flavor = "v2",
    conserve.memory = FALSE,
    return.only.var.genes = TRUE,
    seed.use = 1448145,
    verbose = TRUE,
    ...
){
    set.seed(seed = seed.use)
    vars.to.regress.subset <- vars.to.regress[vars.to.regress %in% colnames(x = object[[]])]
    latent.data <- object[[vars.to.regress.subset]]

    assay = "RNA"
    cell.attr <- slot(object = object, name = 'meta.data')[colnames(object[[assay]]),]
    assay.data <- SCTransform.StdAssay(object = object[[assay]],
                                cell.attr = cell.attr,
                                reference.SCT.model = reference.SCT.model,
                                do.correct.umi = do.correct.umi,
                                ncells = ncells,
                                residual.features = residual.features,
                                variable.features.n = variable.features.n,
                                variable.features.rv.th = variable.features.rv.th,
                                vars.to.regress = vars.to.regress,
                                latent.data = latent.data,
                                do.scale = do.scale,
                                do.center = do.center,
                                clip.range = clip.range,
                                vst.flavor = vst.flavor,
                                conserve.memory = conserve.memory,
                                return.only.var.genes = return.only.var.genes,
                                seed.use = seed.use,
                                verbose = verbose,
                                ...)
    assay.data <- Seurat:::SCTAssay(assay.data, assay.orig = assay)
    
    # Extract all SCT models stored in assay
    sct_models <- slot(object = assay.data, name = "SCTModel.list")
    
    # Update umi.assay field for every SCT model 
    slot(object = assay.data, name = "SCTModel.list") <- lapply(sct_models, function(model) {
        slot(model, name = "umi.assay") <- assay
        model
    })

    object[[new.assay.name]] <- assay.data
    
    return(object)
}


pbmc_data <- Read10X(data.dir = "test_data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc_data)

# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
pbmc <- SCTransform.Seurat(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)


## Checking the result 
# Extract the full matrix
mat <- pbmc@assays$SCT@counts

# Find rows (genes) and columns (cells) with the highest expressions
top_genes <- order(rowSums(mat > 0), decreasing = TRUE)[1:5]
top_cells <- order(colSums(mat > 0), decreasing = TRUE)[1:5]

# Subset the matrix using these top indices
dense_chunk <- mat[top_genes, top_cells]

# Print it as a standard, readable matrix
print(as.matrix(dense_chunk))

# Expected output 
#        TAATGCCTCGTCTC-1 CACCGGGACTTCTA-1 ACGTCGCTCCTGAA-1 CAGTTTACACACGT-1
# RPL15                11                8                6               14
# RPS18                 8                9                2               19
# RPS12                 6                4                4               12
# TMSB4X               47               75               24               62
# RPL10                20               13               17               19
#        GGGAACGAAGCTCA-1
# RPL15                 5
# RPS18                 5
# RPS12                 6
# TMSB4X               68
# RPL10                14