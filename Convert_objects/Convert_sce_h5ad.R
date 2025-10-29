# Aim to convert sce to h5ad

library(SingleCellExperiment)
#library(anndata)
library(reticulate)
library(withr)
reticulate::use_condaenv(condaenv = "convert-objects-env", conda = "/project/ancestryhca/ccohen/mamba_installation/conda/condabin/conda", required = F)
reticulate::import(module = 'anndata')

# read in sce
sce <- readRDS('20231107_09-33_Doublet.dir/RDS_objects.dir/doublet_filtered/SingleCellExp/MSK0785-Ach-Enth_doublet_filtered_SingleCellExp.rds')
sce 

# create an Anndata object
ad <- anndata::AnnData(
    X = t(counts(sce)), # raw matrix, transposed for python
    obs = as.data.frame(colData(sce)),  # cell metadata
    var = as.data.frame(rowData(sce)), # gene metadata
    layers = list(
        counts = t(counts(sce)),
        logcounts = t(logcounts(sce)),
        decontX = t(assays(sce)$decontXcounts),
        soupX = t(assays(sce)$soupX)
    ),
    obsm = list(
        pca = reducedDims(sce)$pca,
        umap = reducedDims(sce)$umap
    )

)

print("Anndata object created")

ad

anndata::write_h5ad(ad, "output.h5ad")

print("object saved")
