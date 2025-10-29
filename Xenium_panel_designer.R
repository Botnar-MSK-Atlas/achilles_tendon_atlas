##################################################
### Prepare datasets for Xenium Panel Designer ###
##################################################


# Datasets used
    # Achilles tendon
    # Quads tendon
    # Bone
    # Developmental
    # Synovium

# Script written using the 10X guide:
# https://www.10xgenomics.com/analysis-guides/creating-single-cell-references-for-xenium-custom-panel-design-from-seurat-or-anndata
# set up

library(tidyverse)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(DropletUtils)

# make a new output folder for each run, with the date & time in the directory name
date <- Sys.Date() %>% str_replace_all("-", "")
time <- format(Sys.time(), "%X") %>% str_replace_all(":", "-") %>%
    str_sub(1, 5)
directory <- paste0(date,"_", time, "_Xenium_panel_designer.dir")
dir.create(directory)


#######################
### Achilles tendon ###
#######################

# read in data
so.achilles <- readRDS("20250210_12-13_Fine_annotation.dir/RDS_objects.dir/Achilles_fine_annotation.rds")

# check counts are integers
# RNA assay
all.equal(GetAssayData(object=so.achilles, assay="RNA", slot="counts")@x, 
          as.integer(GetAssayData(object=so.achilles, assay="RNA", slot="counts")@x)) # TRUE
head(GetAssayData(object=so.achilles, assay="RNA", slot="counts")@x) # Integers
# soupX assay
all.equal(GetAssayData(object=so.achilles, assay="soupX", slot="counts")@x, 
          as.integer(GetAssayData(object=so.achilles, assay="soupX", slot="counts")@x)) # TRUE
head(GetAssayData(object=so.achilles, assay="soupX", slot="counts")@x) # Integers

# get gene names and IDs

mapping <- read.table("Files.dir/mapping.txt", header = TRUE, sep = "\t")

df <- data.frame(hgnc_symbol = rownames(so.achilles))
df <- df %>% left_join(mapping, multiple = "first") %>% 
    # create col with all the ensembl IDs
    mutate(ensembl_gene_id_all = coalesce(ensembl_gene_id, hgnc_symbol)) %>%
    # replace ENS id in hgnc_symbol with NA
    mutate(feature_name = str_replace(hgnc_symbol, "ENSG[:digit:]*", "none"))
df$feature_name <- na_if(df$feature_name, "none")
head(df)

# which metadata do we want
unique(so.achilles$fine_annotation)

# Apply subsampling if required
subsample_rate <- 1

subset <- 1:ncol(so.achilles)
if (subsample_rate < 1) {
    subset <- sample(subset, subsample_rate * length(subset))
}

# create a data matrix and change the rownames to ensemblids
mtx <- GetAssayData(so.achilles, slot = "counts", assay = "soupX")
rownames(mtx) <- df$ensembl_gene_id_all
mtx[1:10, 1:10]

# Create a so with the minimal data
so <- CreateSeuratObject(counts = mtx, 
                         assay = "RNA", # call it RNA even though it might be soupX
                         meta.data = so.achilles[[]])

### Save MEX files with DropletUtils ###
path <- paste0(directory, "/Achilles_tendon")
write10xCounts(
    path,
    GetAssayData(so, assay="RNA", slot="counts")[, subset],
    gene.id = rownames(so),
    gene.symbol = df$feature_name,
    barcodes = colnames(so)[subset],
    type = "sparse",
    version = "3"
)

### Save cell type annotations ###

# Define function
bundleOutputs <- function(out_dir, data, barcodes = colnames(data), 
                          cell_type = "cell_type", subset = 1:length(barcodes)) {
    
    if (require("data.table", quietly = TRUE)) {
        data.table::fwrite(
            data.table::data.table(
                barcode = barcodes,
                annotation = unlist(data[[cell_type]])
            )[subset, ],
            file.path(out_dir, "annotations.csv")
        )
    } else {
        write.table(
            data.frame(
                barcode = barcodes,
                annotation = unlist(data[[cell_type]])
            )[subset, ],
            file.path(out_dir, "annotations.csv"),
            sep = ",", row.names = FALSE
        )
    }
    
    bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))
    
    utils::zip(
        bundle,
        list.files(out_dir, full.names = TRUE),
        zip = "zip"
    )
    
    if (file.info(bundle)$size / 1e6 > 2000) {
        warning("The output file is more than 2G and will need to be subset further.")
    }
}

# Run function
bundleOutputs(out_dir = path, data = so, subset = subset, cell_type = "fine_annotation")

#######################
### Quadriceps tendon ###
#######################

# read in data
so.quads <- readRDS("../20230922_quadriceps/Quads_analysis_RDS/20240523_quads_overall_reannotated.RDS")

# check counts are integers
# RNA assay
all.equal(GetAssayData(object=so.quads, assay="RNA", slot="counts")@x, 
          as.integer(GetAssayData(object=so.quads, assay="RNA", slot="counts")@x)) # TRUE
head(GetAssayData(object=so.quads, assay="RNA", slot="counts")@x) # Integers
# soupX assay
all.equal(GetAssayData(object=so.quads, assay="SoupXcounts", slot="counts")@x, 
          as.integer(GetAssayData(object=so.quads, assay="SoupXcounts", slot="counts")@x)) # TRUE
head(GetAssayData(object=so.quads, assay="SoupXcounts", slot="counts")@x) 
# soupX assay is not integer so have to use RNA assay here

# get gene names and IDs

df <- data.frame(hgnc_symbol = rownames(so.quads))
df <- df %>% left_join(mapping, multiple = "first") %>% 
    # create col with all the ensembl IDs
    mutate(ensembl_gene_id_all = coalesce(ensembl_gene_id, hgnc_symbol)) %>%
    # replace ENS id in hgnc_symbol with NA
    mutate(feature_name = str_replace(hgnc_symbol, "ENSG[:digit:]*", "none"))
df$feature_name <- na_if(df$feature_name, "none")
head(df)

# which metadata do we want
unique(so.quads$cluster_id)

# Apply subsampling if required
subsample_rate <- 1

subset <- 1:ncol(so.quads)
if (subsample_rate < 1) {
    subset <- sample(subset, subsample_rate * length(subset))
}

# create a data matrix and change the rownames to ensemblids
mtx <- GetAssayData(so.quads, slot = "counts", assay = "RNA")
rownames(mtx) <- df$ensembl_gene_id_all
mtx[1:10, 1:10]

# Create a so with the minimal data
so <- CreateSeuratObject(counts = mtx, 
                         assay = "RNA", # call it RNA even though it might be soupX
                         meta.data = so.quads[[]])

### Save MEX files with DropletUtils ###
path <- paste0(directory, "/Quadriceps_tendon")
write10xCounts(
    path,
    GetAssayData(so, assay="RNA", slot="counts")[, subset],
    gene.id = rownames(so),
    gene.symbol = df$feature_name,
    barcodes = colnames(so)[subset],
    type = "sparse",
    version = "3"
)

### Save cell type annotations ###
bundleOutputs(out_dir = path, data = so, subset = subset, cell_type = "cluster_id")

################
### Synovium ###
################

# read in data
so.synovium <- readRDS("../20230302_synovium/20240305_synovium_joints_integrated.RDS")

# check counts are integers
# RNA assay
all.equal(GetAssayData(object=so.synovium, assay="RNA", slot="counts")@x, 
          as.integer(GetAssayData(object=so.synovium, assay="RNA", slot="counts")@x)) # TRUE
head(GetAssayData(object=so.synovium, assay="RNA", slot="counts")@x) # Integers
# soupX assay
all.equal(GetAssayData(object=so.synovium, assay="SoupXcounts", slot="counts")@x, 
          as.integer(GetAssayData(object=so.synovium, assay="SoupXcounts", slot="counts")@x)) # TRUE
head(GetAssayData(object=so.synovium, assay="SoupXcounts", slot="counts")@x) 
# soupX assay is not integer so have to use RNA assay here

# get gene names and IDs

df <- data.frame(hgnc_symbol = rownames(so.synovium))
df <- df %>% left_join(mapping, multiple = "first") %>% 
    # create col with all the ensembl IDs
    mutate(ensembl_gene_id_all = coalesce(ensembl_gene_id, hgnc_symbol)) %>%
    # replace ENS id in hgnc_symbol with NA
    mutate(feature_name = str_replace(hgnc_symbol, "ENSG[:digit:]*", "none"))
df$feature_name <- na_if(df$feature_name, "none")
head(df)

# which metadata do we want
unique(so.synovium$clusterid_general)

# Apply subsampling if required
subsample_rate <- 1

subset <- 1:ncol(so.synovium)
if (subsample_rate < 1) {
    subset <- sample(subset, subsample_rate * length(subset))
}

# create a data matrix and change the rownames to ensemblids
mtx <- GetAssayData(so.synovium, slot = "counts", assay = "RNA")
rownames(mtx) <- df$ensembl_gene_id_all
mtx[1:10, 1:10]

# Create a so with the minimal data
so <- CreateSeuratObject(counts = mtx, 
                         assay = "RNA", # call it RNA even though it might be soupX
                         meta.data = so.synovium[[]])

### Save MEX files with DropletUtils ###
path <- paste0(directory, "/Synovium")
write10xCounts(
    path,
    GetAssayData(so, assay="RNA", slot="counts")[, subset],
    gene.id = rownames(so),
    gene.symbol = df$feature_name,
    barcodes = colnames(so)[subset],
    type = "sparse",
    version = "3"
)

### Save cell type annotations ###
bundleOutputs(out_dir = path, data = so, subset = subset, cell_type = "clusterid_general")


############
### Bone ###
############


# read in the h5ad as sce using zellkonverter
sce.bone <- readH5AD("Files.dir/annotated_cenk_out.h5ad") 
# convert to seurat
so.bone <- as.Seurat(sce.bone, counts = "counts", data = "log1p_norm", project = "Bone")
# change the name of the assay to RNA
so.bone[["RNA"]] <- so.bone[["originalexp"]]
DefaultAssay(so.bone) <- "RNA"

# check counts are integers
# RNA assay
all.equal(GetAssayData(object=so.bone, assay="RNA", slot="counts")@x, 
          as.integer(GetAssayData(object=so.bone, assay="RNA", slot="counts")@x)) # TRUE
head(GetAssayData(object=so.bone, assay="RNA", slot="counts")@x) # Integers
# soupX assay not present

# get gene names and IDs

df <- data.frame(hgnc_symbol = rownames(so.bone))
df <- df %>% left_join(mapping, multiple = "first") %>% 
    # create col with all the ensembl IDs
    mutate(ensembl_gene_id_all = coalesce(ensembl_gene_id, hgnc_symbol)) %>%
    # replace ENS id in hgnc_symbol with NA
    mutate(feature_name = str_replace(hgnc_symbol, "ENSG[:digit:]*", "none"))
df$feature_name <- na_if(df$feature_name, "none")
head(df)

# which metadata do we want
unique(so.bone$Leiden_res0_2_Manual_Cell_Type)

# Apply subsampling if required
subsample_rate <- 1

# create a data matrix and change the rownames to ensemblids
mtx <- GetAssayData(so.bone, slot = "counts", assay = "RNA")
rownames(mtx) <- df$ensembl_gene_id_all
mtx[1:10, 1:10]

# Create a so with the minimal data
so <- CreateSeuratObject(counts = mtx, 
                         assay = "RNA", # call it RNA even though it might be soupX
                         meta.data = so.bone[[]])

### Save MEX files with DropletUtils ###
path <- paste0(directory, "/Bone")
write10xCounts(
    path,
    GetAssayData(so, assay="RNA", slot="counts")[, subset],
    gene.id = rownames(so),
    gene.symbol = df$feature_name,
    barcodes = colnames(so)[subset],
    type = "sparse",
    version = "3"
)

### Save cell type annotations ###
bundleOutputs(out_dir = path, data = so, subset = subset, cell_type = "Leiden_res0_2_Manual_Cell_Type")

#####################
### Developmental ###
#####################


# read in the h5ad as sce using zellkonverter
sce.dev <- readH5AD("../20240116_lattice/developmental/dev_scANVI.h5ad") 

# there are 6 genes with duplicated names
# this appears to update the rownames in both the assays
rownames(sce.dev) <- make.unique(rownames(sce.dev))

# Convert to seurat object
so.dev <- as.Seurat(sce.dev, counts = "counts", 
                    data = "log1p_norm", 
                    project = "Developmental")

# check counts are integers
# RNA assay
all.equal(counts(sce.dev), as.integer(counts(sce.dev))) # TRUE
head(counts(sce.dev))
# soupX assay not present

# get gene names and IDs
# These are present in the rowData
# genes = rownames(sce.dev) or rowData(sce.dev)$Gene
# ensemblIDs = rowData(sce.dev)$ensembl_gene_id

# which metadata do we want
unique(colData(sce.dev)$cell_type)

# Apply subsampling if required
subsample_rate <- 1

# create a data matrix and change the rownames to ensemblids
mtx <- counts(sce.dev)
rownames(mtx) <- rowData(sce.dev)$ensembl_gene_id
mtx[1:10, 1:10]

# change the name of the assay to RNA
so[["RNA"]] <- so[["originalexp"]]
DefaultAssay(so) <- "RNA"

# Create a so with the minimal data
so <- CreateSeuratObject(counts = mtx, 
                         assay = "RNA", # call it RNA even though it might be soupX
                         meta.data = so.dev[[]])

### Save MEX files with DropletUtils ###
path <- paste0(directory, "/Developmental")
write10xCounts(
    path,
    GetAssayData(so, assay="RNA", slot="counts")[, subset],
    gene.id = rownames(so),
    gene.symbol = rownames(sce.dev),
    barcodes = colnames(so)[subset],
    type = "sparse",
    version = "3"
)

### Save cell type annotations ###
bundleOutputs(out_dir = path, data = so, subset = subset, cell_type = "cell_type")
