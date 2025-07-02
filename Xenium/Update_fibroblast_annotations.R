##########################################
#### Achilles annotation update ##########
##########################################

# Update the Achilles annotations
# Use the harmony 0.2 resolution clustering

library(Seurat)

wdir <- "/project/wolf1241/xenium/"

# load the Achilles single cell object
so.achilles <- readRDS(paste0(wdir, "data/Achilles_fine_annotation.rds"))

Idents(so.achilles) <- "soupX_snn_res.0.2"

so.achilles <- RenameIdents(so.achilles, 
                            `0` = "Skeletal muscle cells",
                            `1` = "NEGR1hi fibroblasts",
                            `2` = "COMPhi fibroblasts", 
                            `3` = "Vascular endothelial cells", 
                            `4` = "Skeletal muscle cells", 
                            `5` = "Macrophages", 
                            `6` = "Mural cells", 
                            `7` = "Adipocytes", 
                            `8` = "T cells", 
                            `9` = "Skeletal muscle cells", 
                            `10` = "Lymphatic endothelial cells",
                            `11` = "Satellite cells", 
                            `12` = "Granulocytes", 
                            `13` = "COMPhi fibroblasts",
                            `14` = "B cells", 
                            `15` = "Nervous system cells")

so.achilles$`2fb_annotation` <- Idents(so.achilles)

saveRDS(so.achilles, paste0(wdir, "data/Achilles_fine_annotation.rds"))
