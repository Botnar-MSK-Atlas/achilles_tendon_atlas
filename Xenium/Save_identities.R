# save identities for Xenium explorer

# load objects

so_Enth <- readRDS("20250610_15-13_Xenium_Achilles_Enthesis.dir/RDS_objects.dir/Enthesis_SeuratObject.rds")
so_MB <- readRDS("20250611_11-26_Xenium_Achilles_MB.dir/RDS_objects.dir/MB_SeuratObject.rds")
so_MTJ <- readRDS("20250610_15-27_Xenium_Achilles_MTJ.dir/RDS_objects.dir/MTJ_SeuratObject.rds")
so_muscle <- readRDS("20250610_15-58_Xenium_Achilles_Muscle.dir/RDS_objects.dir/Muscle_SeuratObject.rds")

df <- data.frame(cell_id = colnames(so_Enth),
                 group = so_Enth$predicted.id_2fb)
write.table(df, 
          "20250610_15-13_Xenium_Achilles_Enthesis.dir/Two_fibroblasts_annotation/Enthesis_predicted.id_2fb.csv",
          quote = FALSE, 
          sep = ",",
          col.names = TRUE,
          row.names = FALSE)

df <- data.frame(cell_id = colnames(so_MB),
                 group = so_MB$predicted.id_2fb)
write.table(df, 
            "20250611_11-26_Xenium_Achilles_MB.dir/Two_fibroblasts_annotation/MB_predicted.id_2fb.csv",
            quote = FALSE, 
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)


df <- data.frame(cell_id = colnames(so_MTJ),
                 group = so_MTJ$predicted.id_2fb)
write.table(df, 
            "20250610_15-27_Xenium_Achilles_MTJ.dir/Two_fibroblasts_annotation/MTJ_predicted.id_2fb.csv",
            quote = FALSE, 
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)


df <- data.frame(cell_id = colnames(so_muscle),
                 group = so_muscle$predicted.id_2fb)
write.table(df, 
            "20250610_15-58_Xenium_Achilles_Muscle.dir/Two_fibroblasts_annotation/Muscle_predicted.id_2fb.csv",
            quote = FALSE, 
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)

