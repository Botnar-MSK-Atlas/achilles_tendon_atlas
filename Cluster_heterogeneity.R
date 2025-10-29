# Cluster heterogeneity using ROGUE

# Score cluster heterogeneity using ROGUE. 
# https://github.com/PaulingLiu/ROGUE
# https://www.nature.com/articles/s41467-020-16904-3


library(tidyverse)
library(Seurat)
library(ROGUE)

# make a new output folder for each run, with the date & time in the directory name
date <- Sys.Date() %>% str_replace_all("-", "")
time <- format(Sys.time(), "%X") %>% str_replace_all(":", "-") %>%
    str_sub(1, 5)
directory <- paste0(date,"_", time, "_Cluster_heterogeneity.dir")
dir.create(directory, showWarnings = FALSE)

# Read in the Achilles subset fibroblast object and the full Achilles object. 

so.fibroblast <- readRDS("20241008_14-55_Subset_fibroblasts.dir/RDS_objects.dir/Achilles_fibroblast_subset.rds")
so.full <- readRDS("20240917_11-29_Annotation_update.dir/RDS_objects.dir/Achilles_integrated_annotated.rds")

# subset the fibroblasts from the full object
Idents(so.full) <- so.full$cell_annotation_scVI_0.2
so.full.subset <- subset(so.full, idents = c("NEGR1+ fibroblasts", "ITGA10+ fibroblasts"))

# Make a function to run the ROGUE on a seurat object

run_rogue <- function(so, name){
    
    # Extract the expression matrix & metadata
    expr <- as.matrix(GetAssayData(so, assay = "soupX", slot = "data"))
    meta <- so[[]]
    
    # Filter out low quality genes and low quality cells
    expr.filter <- matr.filter(expr, min.cells = 10, min.genes = 10)
    
    # Apply the expression entropy model & plot
    
    ent.res <- SE_fun(expr.filter)
    SEplot(ent.res)
    ggsave(ggsave(paste0(directory, "/SE_plot_", name, ".png")))
    
    # Calculate ROGUE score across all cells
    
    rogue.value <- CalculateRogue(ent.res, platform = "UMI")
    
    # Calculate and plot ROGUE scores for all calculated resolutions
    
    resolutionList <- grep("leiden_X_scVI_snn_res", colnames(meta), value = TRUE)
    
    rogue.res.list <- list()
    colmeans.list <- list()
    dir.create(paste0(directory, "/ROGUE_scores_", name))
    dir.create(paste0(directory, "/Boxplots_", name))
    for (resolution in resolutionList){
        print(resolution)
        rogue.res.list[[resolution]] <- rogue(expr.filter, labels = meta[[resolution]], 
                                              samples = meta$patient, platform = "UMI", span = 0.6)
        write.table(rogue.res.list[[resolution]], 
                    paste0(directory, "/ROGUE_scores_", name, "/", resolution, ".txt"),
                    quote = FALSE)
        colmeans.list[[resolution]] <- colMeans(rogue.res.list[[resolution]], na.rm = T)
        rogue.boxplot(rogue.res.list[[resolution]])+
            geom_hline(yintercept = rogue.value, linetype = "dashed")+
            ggtitle(resolution)
        ggsave(paste0(directory, "/Boxplots_", name, "/", resolution, ".png"))
    }
    
    # Calculate and plot the mean ROGUE scores per cluster per sample
    
    df <- plyr::ldply(colmeans.list, rbind) %>% 
        rename(resolution = .id)
    
    df <- df %>% as_tibble() %>%
        pivot_longer(cols = 2:ncol(df), names_to = "cluster", values_to = "score")
    
    # make cluster a factor
    df$cluster <- as.numeric(df$cluster)
    df$cluster <- factor(df$cluster, levels = sort(unique(df$cluster)))
    
    ggplot(df, aes(x = cluster, y = score, group = resolution))+
        geom_line(na.rm = T, aes(colour = resolution))+
        geom_point(na.rm = T, aes(colour = resolution))+
        theme_classic()
    ggsave(paste0(directory, "/Mean_score_per_cluster_", name, ".png"), width = 10, height = 5)
    
    # Calculate and plot the mean ROGUE scores per resolution
    
    df <- plyr::ldply(colmeans.list, rbind) %>% 
        rename(resolution = .id)
    # move resolution to rownames
    rownames(df) <- df$resolution
    df$resolution <- NULL
    
    # how many colums are populated per row?
    my_fun <- function(x){
        sum(!is.na(x))
    }
    row_means <- data.frame("resolution" = str_remove(resolutionList, "leiden_X_scVI_snn_res."), 
                            "mean_score" = rowMeans(df, na.rm = T), 
                            "n_clusters" = apply(df, 1, my_fun))
    
    ggplot(row_means, aes(x = mean_score, y = resolution))+
        geom_point(aes(size = n_clusters))+
        theme_classic()
    ggsave(paste0(directory, "/Mean_score_per_resolution_", name, ".png"))
    
}

run_rogue(so.full, "Achilles")
run_rogue(so.fibroblast, "Fibroblast_subset_reintegrated")
run_rogue(so.full.subset, "Achilles_fibroblast_subset")


