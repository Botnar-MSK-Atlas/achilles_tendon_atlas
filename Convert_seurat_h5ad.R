# Aim to convert Seurat to h5ad

library(Seurat)
#library(anndata)
library(reticulate)
library(withr)
library(tidyverse)
library(Matrix)
reticulate::use_condaenv(condaenv = "convert-objects-env", conda = "/project/ancestryhca/ccohen/mamba_installation/conda/condabin/conda", required = F)
reticulate::import(module = 'anndata')

# make a new output folder for each run, with the date & time in the directory name
date <- Sys.Date() %>% str_replace_all("-", "")
time <- format(Sys.time(), "%X") %>% 
    str_replace_all(":", "-") %>%
    str_sub(1,5)
directory <- paste0(date,"_", time, "_convert-objects.dir")
dir.create(directory)

# read in Seurat object
so <- readRDS("20241008_12-04_Integration.dir/RDS_objects.dir/Achilles_harmony_SeuratObject.rds")
so 

DefaultAssay(so) <- "RNA"
# create an Anndata object
ad <- anndata::AnnData(
    X = t(GetAssayData(so, slot = "counts")), # raw matrix, transposed for python
    obs = as.data.frame(so[[]]),  # cell metadata
    var = as.data.frame(rownames(so)), # gene metadata
    layers = list(
        counts = t(GetAssayData(so, slot = "counts")),
        logcounts = t(GetAssayData(so, slot = "data")),
        #decontX = t(so@assays$decontXcounts$counts),
        soupX = t(so@assays$soupX@counts)
    ),
    obsm = list(
        pca = Embeddings(so, reduction = "lognorm.pca")#,
        #umap = Embeddings(so, reduction = "harmony.umap")
    )
    
)

print("Anndata object created")

ad

anndata::write_h5ad(ad, paste0(directory, "/Achilles_integrated_annotated.h5ad"))

print("object saved")
