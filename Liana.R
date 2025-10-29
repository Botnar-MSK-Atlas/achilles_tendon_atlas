#############
### LIANA ###
#############

# Run Liana on Achilles data 

library(tidyverse)
library(Seurat)
library(liana)

# make a new output folder for each run, with the date & time in the directory name
date <- Sys.Date() %>% str_replace_all("-", "")
time <- format(Sys.time(), "%X") %>% str_replace_all(":", "-") %>%
    str_sub(1, 5)
directory <- paste0(date,"_", time, "_Liana.dir")
dir.create(directory, showWarnings = FALSE)

# read in the data
so.achilles <- readRDS("20250228_12-50_Fine_annotation.dir/RDS_objects.dir/Achilles_fine_annotation.rds")
# downsample for testing
# so.achilles <- subset(so.achilles, downsample = 100)

# Make a mid-level annotation with fine annotation of fibroblasts and broad annotation of other cell types

metadata <- so.achilles[[]] %>% select(fine_annotation, broad_annotation)
df <- data.frame(fine_annotation = unique(metadata$fine_annotation))
df$fine_annotation_short <- df$fine_annotation %>% 
    str_replace("Nervous system cells", "Schwann cells") %>% 
    str_replace("NEGR1hi VCANhi fibroblasts", "NEGR1.VCAN.fb") %>% 
    str_replace("NEGR1hi ITGA6hi fibroblasts", "NEGR1.ITGA6.fb") %>% 
    str_replace("COMPhi THBS4hi fibroblasts"  , "COMP.THBS4.fb") %>% 
    str_replace("PRG4hi fibroblasts", "PRG4.fb") %>% 
    str_replace("NEGR1hi COL15A1hi fibroblasts", "NEGR1.COL15A1.fb") %>% 
    str_replace("COMPhi MMP3hi fibroblasts" , "COMP.MMP3.fb") %>% 
    str_replace("Lymphatic endothelial cells", "LEC") %>% 
    str_replace("Slow-twitch skeletal muscle cells", "Slow muscle") %>%
    str_replace("Fast-twitch skeletal muscle cells", "Fast muscle") %>%
    str_replace("Transitional skeletal muscle cells", "Transitional muscle")

df$fibroblast_annotation <- df$fine_annotation_short %>% 
    str_replace("MERTKhi LYVE1hi macrophages", "Macrophages") %>%
    str_replace("MERKhi LYVE1lo macrophages", "Macrophages") %>%
    str_replace("MERTKlo PTPRGhi macrophages", "Macrophages") %>%
    str_replace("CLEC10Ahi DCs", "Macrophages") %>%
    str_replace("Monocytes", "Macrophages") %>%
    str_replace("CLEC9Ahi DCs", "Macrophages") %>%
    str_replace("Venular VEC", "VEC") %>%
    str_replace("Arteriolar VEC", "VEC") %>%
    str_replace("vSMC", "Mural cells") %>%
    str_replace("Pericytes", "Mural cells") %>%
    str_replace("Slow muscle", "Skeletal muscle cells") %>%
    str_replace("Fast muscle", "Skeletal muscle cells") %>%
    str_replace("Transitional muscle", "Skeletal muscle cells") %>%
    str_replace("Plasma cells", "B cells") 
    
metadata <- left_join(metadata, df)
so.achilles <- AddMetaData(so.achilles, metadata$fibroblast_annotation, "fibroblast_annotation")
so.achilles <- AddMetaData(so.achilles, metadata$fine_annotation_short, "fine_annotation_short")

so.achilles[[]] %>% select(fine_annotation, fibroblast_annotation, fine_annotation_short) %>% head(20)

# Run Liana using 7 methods with fibroblast annotation
Idents(so.achilles) <- so.achilles$fibroblast_annotation

liana_test <- liana_wrap(so.achilles)
liana_test <- liana_test %>%
    liana_aggregate()

# save results

write.table(liana_test, paste0(directory, "/Liana_aggregate_results_fibroblast_annotation.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


# Run Liana using 7 methods with fine annotation
Idents(so.achilles) <- so.achilles$fine_annotation_short

liana_test <- liana_wrap(so.achilles)
liana_test <- liana_test %>%
    liana_aggregate()

# save results

write.table(liana_test, paste0(directory, "/Liana_aggregate_results_fine_annotation.txt"), sep = "\t", row.names = FALSE, quote = FALSE)





