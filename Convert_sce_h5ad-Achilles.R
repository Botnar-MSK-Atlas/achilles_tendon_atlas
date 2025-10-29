# Convert sce to h5ad for Achilles dataset

library(SingleCellExperiment)
library(reticulate)
library(withr)
library(tidyverse)
reticulate::use_condaenv(condaenv = "convert-objects-env", conda = "/project/ancestryhca/ccohen/mamba_installation/conda/condabin/conda", required = F)
reticulate::import(module = 'anndata')

# make a new output folder for each run, with the date & time in the directory name
date <- Sys.Date() %>% str_replace_all("-", "")
time <- format(Sys.time(), "%X") %>% 
    str_replace_all(":", "-") %>%
    str_sub(1,5)
directory <- paste0(date,"_", time, "_convert-objects.dir")
dir.create(directory, showWarnings = FALSE)

# read in sce
#make a list of sample files
file_list <- list.files("20240827_10-27_Doublet.dir/RDS_objects.dir/doublet_filtered/SingleCellExp/", pattern=".rds", recursive = TRUE)
names(file_list) <- str_replace(file_list, "_doublet_filtered_SingleCellExp.rds", "") 

#Make a new list for the sce objects
sce <- list()

# generate a list of sce objects
for (i in 1:length(file_list)){
    sce[[i]] <- readRDS(paste0("20240827_10-27_Doublet.dir/RDS_objects.dir/doublet_filtered/SingleCellExp/",file_list[i]))
}
sce

# function to create an anndata object
create_anndata <- function(sce){
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
    return(ad)
}

# convert sce to anndata
h5ad <- list()

for (i in 1:length(sce)){
    h5ad[[i]] <- create_anndata(sce[[i]])
}

h5ad


print("Anndata objects created")

# save objects
dir.create(paste0(directory, "/h5ad_objects.dir/"))
for (i in 1:length(h5ad)){
    anndata::write_h5ad(h5ad[[i]], paste0(directory, "/h5ad_objects.dir/", names(file_list)[[i]], ".h5ad"))
}
print("objects saved")
