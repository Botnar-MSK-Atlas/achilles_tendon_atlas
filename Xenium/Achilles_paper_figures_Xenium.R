##############################
# Xenium Figures for Achilles paper #
##############################

#--------
# set up
#--------

library(Seurat)
library(tidyverse)
library(scCustomize)
library(readr)
library(pheatmap)
library(matrixStats)
library(spdep)
library(geojsonR)
library(clustree)
library(yaml)
library(cowplot)
library(RColorBrewer)
library(SCpubr)

wdir <- "/project/wolf1241/xenium/"
set.seed(123)

#------------------------
# make output directories
#------------------------

directory <- "Achilles_paper_figures.dir"
dir.create(directory)
dir.create(paste0(directory, "/Figure_1/"))
dir.create(paste0(directory, "/Figure_2/"))
dir.create(paste0(directory, "/Figure_3/"))
dir.create(paste0(directory, "/Supplementary"))
dir.create(paste0(directory, "/Supplementary/Figure_S4/"))
dir.create(paste0(directory, "/Supplementary/Figure_S5/"))
dir.create(paste0(directory, "/Supplementary/Figure_S6/"))
dir.create(paste0(directory, "/Supplementary/Figure_S7/"))
dir.create(paste0(directory, "/Supplementary/Figure_S8/"))
#------------------------
# set the colours
#------------------------



broad.colours <- c("COMPhi MMP3hi fibroblasts" =  "#F9BA73",
                   "COMPhi THBS4hi fibroblasts" = "#44BED0",
                   "NEGR1hi COL15A1hi fibroblasts" = "#41A021" ,
                   "NEGR1hi VCANhi fibroblasts" = "#FC6400",
                   "NEGR1hi ITGA6hi fibroblasts" = "#C19C93",
                   "Schwann cells" = "#458BFF",
                   "Immune cells" = "#44BED0",
                   "Endothelial cells"  ="#FFABC3" ,
                   "Adipocytes" = "#DADB89" ,
                   "Skeletal muscle cells" = "#9267BF",
                   "Mural cells" = "#7CEC37",
                   "Chondrocytes" = "#C7C7C7", 
                   "PRG4hi fibroblasts" = "#7F7F7F"
)  


achilles.colours<- c(# Broad groups
  "Fibroblasts" = "#FC6400", 
  "Skeletal muscle cells" = "#9267BF", 
  "Satellite cells" = "#C4B0D6",
  "Mural cells" = "#7CEC37",
  # Immune
  "Macrophages" = "#44BED0", 
  "B cells" = "#88564A", 
  "T cells" = "#41A021", 
  "Plasma cells" = "#01EBD6", 
  "Granulocytes" ="#2D5780",
  # Stromal
  "Vascular endothelial cells" = "#FFABC3",
  "Lymphatic endothelial cells" = "#DD51AD",
  "Adipocytes" = "#DADB89", 
  "Mural" = "#7CEC37",
  "Schwann cells" = "#458BFF", 
  #"Nervous system cells" = "#F67E00", 
  "Ambient" = "grey")


# Fibroblast colours 
fibroblast.colours <- c("COMPhi MMP3hi fibroblasts" =  "#F9BA73",
                        "COMPhi THBS4hi fibroblasts" = "#44BED0",
                        "NEGR1hi COL15A1hi fibroblasts" = "#41A021" ,
                        "NEGR1hi VCANhi fibroblasts" = "#FC6400",
                        "NEGR1hi ITGA6hi fibroblasts" = "#C19C93",
                        "PRG4hi fibroblasts"= "#67001F",
                        "Chondrocytes" = "#C7C7C7"
)


# Muscle colours

muscle.colours <- c("Slow-twitch skeletal muscle cells" = "#9267BF",
                    "Transitional skeletal muscle cells" = "#44BED0", 
                    "Satellite cells" = "#CD2321", 
                    "Fast-twitch skeletal muscle cells" = "#BBBD00")


# Immune colours

immune.colours <- c("B cells" = "#88564A", 
                    "Monocytes" = "#21EFA5",
                    "T cells" = "#FA2CF6", 
                    "CLEC10Ahi DCs" = "#E8EF21",
                           "Granulocytes" = "#CD2321",
                    "MERTKhi LYVE1hi macrophages" = "#43BCE7",
                    "NK cells" = "#F9BA73"
                           )


# stromal colours

stromal.colours <- c("Lymphatic endothelial cells" = "#FE9A30", 
                     "Adipocytes" = "#D2C800",
                     "Arteriolar VEC" = "#F2B6D2",
                     "Venular VEC" =  "#DC77C3",
                     "vSMC" = "#9FDF86",
                     "Pericytes" = "#41A021",
                     "Schwann cells" = "#458BFF"
)



ma.cols <-  c(Enthesis = "#d37504", MB = "#959295", MTJ = "#359c72", Muscle = "#34688d")

seaborn.colours <- c("#3378B6", "#B1C7E9",
                     "#F67E00", "#F9BA73",
                     "#41A021", "#9FDF86", 
                     "#CD2321", "#F89795", 
                     "#9267BF", "#C4B0D6",
                     "#88564A", "#C19C93",
                     "#DC77C3", "#F2B6D2",
                     "#7F7F7F", "#C7C7C7", 
                     "#BBBD00", "#DADB89", 
                     "#44BED0", "#A5DAE6")

#------------------------
# Read in the objects
#------------------------

so.harmony <- readRDS(paste0(wdir, "20250708_10-26_Visualisation_.dir/RDS_objects.dir/Achilles_integrated_annotated.rds"))
so.achilles <- readRDS(paste0(wdir, "data/Achilles_fine_annotation.rds"))

so.fb <- readRDS(paste0(wdir, "20250702_12-14_Fibroblast_annotation.dir/RDS_objects.dir/Achilles_fibroblast_subset.rds"))
so.im <- readRDS(paste0(wdir, "20250702_13-21_Immune_annotation_.dir/RDS_objects.dir/Achilles_immune_subset.rds"))
so.str <- readRDS(paste0(wdir, "20250723_11-23_Stromal_annotation_.dir/RDS_objects.dir/Achilles_stromal_subset.rds"))
so.musc <- readRDS(paste0(wdir, "20250723_11-47_Muscle_annotation_.dir/RDS_objects.dir/Achilles_muscle_subset.rds"))

# rename nervous system cells to schwann cells
# rename chondrocytes to COMPhi MMP3hi fibroblasts

Idents(so.achilles) <- "broad_annotation"
so.achilles <- RenameIdents(so.achilles, 
                       'Nervous system cells' = "Schwann cells")
so.achilles$broad_annotation <- Idents(so.achilles)


Idents(so.harmony) <- "predicted.id_integrated_broad"
so.harmony <- RenameIdents(so.harmony, 
                           `Nervous system cells` = "Schwann cells")
so.harmony$predicted.id_integrated_broad <- Idents(so.harmony)

Idents(so.str) <- "predicted.id_stromal_fine"
so.str <- RenameIdents(so.str, 
                           'Nervous system cells' = "Schwann cells")
so.str$predicted.id_stromal_fine <- Idents(so.str)

Idents(so.fb) <- "predicted.id_fb_fine"
so.fb <- RenameIdents(so.fb, 
                      `Chondrocytes` = "COMPhi MMP3hi fibroblasts")
so.fb$predicted.id_fb_fine <- Idents(so.fb)

# Update the fine annotations in the integrated object

df <- so.harmony[[]]
df.fb <- so.fb[[]] %>% select(cell_id, predicted.id_fb_fine)
df.im <- so.im[[]] %>% select(cell_id, predicted.id_immune_fine)
df.musc <- so.musc[[]] %>% select(cell_id, predicted.id_muscle_fine)
df.str <- so.str[[]] %>% select(cell_id, predicted.id_stromal_fine)
df <- left_join(df, df.fb)
df <- left_join(df, df.im)
df <- left_join(df, df.musc)
df <- left_join(df, df.str)
df <- df %>% select(broad_annotation, predicted.id_broad, predicted.id_fb_fine, predicted.id_immune_fine, predicted.id_muscle_fine, predicted.id_stromal_fine)
head(df)

# make a new column of fine annotation by filling in the values from the other annotations
df <- df %>% mutate(fine_annotation = coalesce(predicted.id_fb_fine, 
                                               predicted.id_immune_fine, 
                                               predicted.id_muscle_fine, 
                                               predicted.id_stromal_fine))

# make a new column of fibroblast annotaiton by using fb fine annotation and everything else broad annotation

df <- df %>% mutate(fibroblast_annotation = coalesce(predicted.id_fb_fine, 
                                                     broad_annotation))

# add these columns to the object
so.harmony$fine_annotation <- df$fine_annotation
so.harmony$fibroblast_annotation <- df$fibroblast_annotation

#-----------------------------------------------------
# Figure 1E: Export the annotations for Xenium Explorer
#-----------------------------------------------------


# Broad annotation (same as single nuc)
df <- data.frame(cell_id = colnames(so.harmony),
                 group = so.harmony$predicted.id_broad)
write.table(df, 
            paste0(directory, "/Figure_1/Broad_annotation.csv"),
            quote = FALSE, 
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)

# Fine annotation
df <- data.frame(cell_id = colnames(so.harmony),
                 group = so.harmony$fine_annotation)
write.table(df, 
            paste0(directory, "/Figure_1/Fine_annotation.csv"),
            quote = FALSE, 
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)

# Fibroblast annotation
df <- data.frame(cell_id = colnames(so.harmony),
                 group = so.harmony$fibroblast_annotation)
write.table(df, 
            paste0(directory, "/Figure_1/Fibroblast_annotation.csv"),
            quote = FALSE, 
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)


#----------------------------------------------------------
# Figure 2F: Distributions of COMP and NEGR1 across the tendon
#----------------------------------------------------------

# Subset the so to the MB
so.mb <- subset(so.harmony, microanatomical_site == "MB")
so.mb[[]] <- so.mb[[]] %>% mutate (coords_sum = x_centroid + y_centroid)

# subset to cells in the slice
so.mb$cells_slice <- so.mb$coords_sum > 3900 & so.mb$coords_sum < 4100
so.slice <- subset(so.mb, cells_slice)
# ImageDimPlot(so.slice, group.by = "broad_annotation")

df$COMP <- GetAssayData(so.slice, assay = "SCT", layer = "scale.data")["COMP",]
df$NEGR1 <- GetAssayData(so.slice, assay = "SCT", layer = "scale.data")["NEGR1",]

df <- df %>% filter (broad_annotation == 'COMPhi MMP3hi fibroblasts' |
                     broad_annotation == "NEGR1hi VCANhi fibroblasts" |
                     broad_annotation == "NEGR1hi ITGA6hi fibroblasts" |  
                     #broad_annotation == "COMPhi THBS4hi fibroblasts" |
                     broad_annotation == "NEGR1hi COL15A1hi fibroblasts") 

p1 <- ggplot(df, aes(x = x_centroid, y = COMP))+
  geom_point(aes(colour = broad_annotation))+
  #geom_segment( aes(x=x_centroid, xend=x_centroid, y=0, yend=celltype), color="grey") +
  scale_colour_manual(values = broad.colours)+
  theme_classic()+
  geom_smooth(colour = "black")+
  ylab("scaled expression")+
  ggtitle("COMP")+
  theme(legend.position="none")


p2 <- ggplot(df, aes(x = x_centroid, y = NEGR1))+
  geom_point(aes(colour = broad_annotation))+
  #geom_segment( aes(x=x_centroid, xend=x_centroid, y=0, yend=celltype), color="grey") +
  scale_colour_manual(values = broad.colours)+
  theme_classic()+
  geom_smooth(colour = "black")+
  ylab("scaled expression")+
  ggtitle("NEGR1")+
  theme(legend.position="none")


plot_grid(p1, p2, ncol = 1)

ggsave(paste0(directory, "/Figure_2/COMP_NEGR1_across_slice.png"), 
       width = 8, height = 6, bg = "white")

#--------------------------------------------------
# Figure S4A: UMAP of predicted broad annotations
#---------------------------------------------------

do_DimPlot(so.harmony, reduction = "harmony.umap", group.by = "predicted.id_integrated_broad",
           colors.use = achilles.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right"
)


ggsave(paste0(directory, "/Supplementary/Figure_S4/UMAP_predicted_broad.png"), width = 10, height = 7, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S4/UMAP_predicted_broad.svg"), width = 10, height = 7, bg = "white")


#--------------------------------------------------
# Figure S4B: Cell proportions in snRNAseq vs Xenium
#---------------------------------------------------

# what are the cell numbers for broad annotations ?

df <- so.harmony[[]] %>% select("predicted.id_integrated_broad", "fine_annotation")

df.broad <- df %>% 
  group_by(predicted.id_integrated_broad) %>%  
  summarise(xenium=n()) %>% 
  mutate(xenium_percent=(xenium/sum(xenium))*100) %>% 
  arrange(desc(xenium))

df.broad$predicted.id_integrated_broad <- factor(df.broad$predicted.id_integrated_broad, 
                                                 levels = df.broad$predicted.id_integrated_broad)

df.sn <- so.achilles[[]] %>% select(broad_annotation)

df.sn <- df.sn %>% 
  group_by(broad_annotation) %>%  
  summarise(snRNAseq=n()) %>% 
  mutate(snRNAseq_percent=(snRNAseq/sum(snRNAseq))*100) %>% 
  arrange(desc(snRNAseq))
df.sn$broad_annotation <- factor(df.sn$broad_annotation, levels = df.sn$broad_annotation)

df <- left_join(df.sn, df.broad, by = c("broad_annotation" = "predicted.id_integrated_broad"))

df.proportion <- df %>% select(broad_annotation, snRNAseq_percent, xenium_percent)
df.proportion$broad_annotation <- factor(df.proportion$broad_annotation, levels = df.proportion$broad_annotation)
df.proportion <- df.proportion %>% pivot_longer(2:3, names_to = "technology", values_to = "percent_cells")

ggplot(df.proportion, aes(y=percent_cells, x=broad_annotation, fill = technology))+
  geom_bar(position = "dodge", stat = "identity")+
  theme_classic()+
  scale_fill_manual(values = c("snRNAseq_percent" = "#d95f02", "xenium_percent" = "#7570b3"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank())

ggsave((paste0(directory, "/Supplementary/Figure_S4/Cell_percent.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S4/Cell_percent.svg")), 
       width = 8, height = 5)

#--------------------------------------------------
# Figure S4D: Cell proportions across microanatomy
#--------------------------------------------------- 


df <- so.harmony[[]] %>% 
  select("predicted.id_integrated_broad", "microanatomical_site")%>% 
  group_by(predicted.id_integrated_broad, microanatomical_site) %>%  
  summarise(cell_number=n()) %>% 
  mutate(percent_cells=(cell_number/sum(cell_number))*100)


df$predicted.id_integrated_broad <- factor(df$predicted.id_integrated_broad, 
                              levels = df.broad$predicted.id_integrated_broad)

ggplot(df, aes(y=percent_cells, x=predicted.id_integrated_broad, fill = microanatomical_site))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_manual(values = ma.cols)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x = "", y = "Proportion")+
  labs(fill = "Microanatomical site")

ggsave((paste0(directory, "/Supplementary/Figure_S4/Xenium_broad_annotation_microanatomy_percent.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S4/Xenium_broad_annotation_microanatomy_percent.svg")), 
       width = 8, height = 5)

#--------------------------------------------------
# Figure S5D: Fibroblast fine annotation UMAP
#---------------------------------------------------


do_DimPlot(so.fb, reduction = "umap", group.by = "predicted.id_fb_fine",
           colors.use = broad.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right"
)

ggsave(paste0(directory, "/Supplementary/Figure_S5/UMAP_predicted_fibroblast.png"), 
       width = 8, height = 4, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S5/UMAP_predicted_fibroblast.svg"), 
       width = 8, height = 4, bg = "white")


#--------------------------------------------------
# Figure S5I: Fibroblast cell proportions in snRNAseq and Xenium
#---------------------------------------------------

Idents(so.harmony) <- so.harmony$predicted.id_integrated_broad
so.fb <- subset(so.harmony, idents = "Fibroblasts")
df <- so.fb[[]] %>% select("broad_annotation", "fine_annotation")

df.fine <- df %>% 
  group_by(fine_annotation) %>%  
  summarise(xenium=n()) %>% 
  mutate(xenium_percent=(xenium/sum(xenium))*100) %>% 
  arrange(desc(xenium))

df.fine$fine_annotation <- factor(df.fine$fine_annotation, 
                                       levels = df.fine$fine_annotation)

Idents(so.achilles) <- so.achilles$broad_annotation
so.sn.fb <- subset(so.achilles, idents = "Fibroblasts")
df.sn <- so.sn.fb[[]] %>% select(fine_annotation)

df.sn <- df.sn %>% 
  group_by(fine_annotation) %>%  
  summarise(snRNAseq=n()) %>% 
  mutate(snRNAseq_percent=(snRNAseq/sum(snRNAseq))*100) %>% 
  arrange(desc(snRNAseq))
df.sn$fine_annotation <- factor(df.sn$fine_annotation, 
                                levels = df.sn$fine_annotation)

df <- full_join(df.fine, df.sn)
df.proportion <- df %>% select(fine_annotation, snRNAseq_percent, xenium_percent)
df.proportion <- df.proportion %>% pivot_longer(2:3, names_to = "technology", values_to = "percent_cells")

ggplot(df.proportion, aes(y=percent_cells, x=fine_annotation, fill = technology))+
  geom_bar(position = "dodge", stat = "identity")+
  theme_classic()+
  scale_fill_manual(values = c("snRNAseq_percent" = "#d95f02", "xenium_percent" = "#7570b3"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank())

ggsave((paste0(directory, "/Supplementary/Figure_S5/Fibroblast_cell_percent.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S5/Fibroblast_cell_percent.svg")), 
       width = 8, height = 5)

# Cell proportions for broad annotation split by microanatomy


df <- so.fb[[]] %>% 
  select("predicted.id_fb_fine", "microanatomical_site")%>% 
  group_by(predicted.id_fb_fine, microanatomical_site) %>%  
  summarise(cell_number=n()) %>% 
  mutate(percent_cells=(cell_number/sum(cell_number))*100)


df$predicted.id_fb_fine <- factor(df$predicted.id_fb_fine, 
                                  levels = df.fine$predicted.id_fb_fine)

ggplot(df, aes(y=cell_number, x=predicted.id_fb_fine, fill = microanatomical_site))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_manual(values = ma.cols)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank())

ggsave((paste0(directory, "/Supplementary/Figure_S/Xenium_fibroblast_annotation_microanatomy.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S5/Xenium_fibroblast_annotation_microanatomy.svg")), 
       width = 8, height = 5)

ggplot(df, aes(y=percent_cells, x=predicted.id_fb_fine, fill = microanatomical_site))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_manual(values = ma.cols)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x = "", y = "Proportion")+
  labs(fill = "Microanatomical site")

ggsave((paste0(directory, "/Supplementary/Figure_S5/Xenium_fibroblast_microanatomy_percent.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S5/Xenium_fibroblast_microanatomy_percent.svg")), 
       width = 8, height = 5)



#--------------------------------------------------
# Figure S6D: Stromal fine annotation
#---------------------------------------------------


do_DimPlot(so.str, reduction = "umap", group.by = "predicted.id_stromal_fine",
           colors.use = stromal.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right"
)

ggsave(paste0(directory, "/Supplementary/Figure_S7/UMAP_predicted_stromal.png"), 
       width = 8, height = 4, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S7/UMAP_predicted_stromal.svg"), 
       width = 8, height = 4, bg = "white")



#--------------------------------------------------
# Figure S7D: Immune fine annotation
#---------------------------------------------------

do_DimPlot(so.im, reduction = "umap", group.by = "predicted.id_immune_fine",
           colors.use = immune.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right"
)

ggsave(paste0(directory, "/Supplementary/Figure_S8/UMAP_predicted_immune.png"), 
       width = 8, height = 4, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S8/UMAP_predicted_immune.svg"), 
       width = 8, height = 4, bg = "white")

#--------------------------------------------------
# Figure S8D: Muscle fine annotation
#---------------------------------------------------


do_DimPlot(so.musc, reduction = "umap", group.by = "predicted.id_muscle_fine",
           colors.use = muscle.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right"
)

ggsave(paste0(directory, "/Supplementary/Figure_S6/UMAP_predicted_muscle.png"), 
       width = 8, height = 4, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S6/UMAP_predicted_muscle.svg"), width = 10, height = 7, bg = "white")


#-------------------------------------------------------
# Record session Info
#-------------------------------------------------------


sink(paste0(directory, "/sessionInfo.txt"))
sessionInfo()
sink()


