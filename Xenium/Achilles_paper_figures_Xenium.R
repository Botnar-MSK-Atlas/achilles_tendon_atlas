##############################
# Figures for Achilles paper #
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
dir.create(paste0(directory, "/Supplementary"))
dir.create(paste0(directory, "/Supplementary/Figure_S4/"))
dir.create(paste0(directory, "/Supplementary/Figure_S5/"))

#------------------------
# set the colours
#------------------------

broad.colours <- c("COMPhi MMP3hi fibroblasts" =  "#F9BA73",
                   "COMPhi THBS4hi fibroblasts" = "#44BED0",
                   "NEGR1hi COL15A1hi fibroblasts" = "#41A021" ,
                   "NEGR1hi VCANhi fibroblasts" = "#CD2321",
                   "NEGR1hi ITGA6hi fibroblasts" = "#C19C93",
                   "Nervous system cells" = "#F67E00",
                   "Immune cells" = "#A5DAE6",
                   "Endothelial cells"  ="#DC77C3" ,
                   "Adipocytes" = "#DADB89" ,
                   "Skeletal muscle cells" = "#9267BF",
                   "Mural cells" = "#9FDF86",
                   "Chondrocytes" = "#C7C7C7", 
                   "PRG4hi fibroblasts" = "#7F7F7F"
)  

achilles.colours <- c(# Broad groups
  "Fibroblasts" = "#B2182B", 
  "Skeletal muscle cells" = "#88419D", 
  "Satellite cells" = "#8C96C6",
  "Mural cells" = "#7FBC41",
  # Immune
  "Macrophages" = "#35978F", 
  "B cells" = "#543005", 
  "T cells" = "#F6E8C3", 
  "Plasma cells" = "#01EBD6", 
  "Granulocytes" = "#80CDC1", 
  # Stromal
  "Vascular endothelial cells" = "#F1B6DA",
  "Lymphatic endothelial cells" = "#C51B7D",
  "Adipocytes" = "#E6F598", 
  "Mural" = "#B8E186", 
  "Nervous system cells" = "#276419", 
  "Ambient" = "grey")

# Muscle colours
muscle.colours <- brewer.pal(5, "BuPu")
names(muscle.colours) <- c("unused", 
                           "Transitional skeletal muscle cells", 
                           "Satellite cells", 
                           "Fast-twitch skeletal muscle cells", 
                           "Slow-twitch skeletal muscle cells")


ma.cols <-  c(Enthesis = "#d37504", MB = "#959295", MTJ = "#359c72", Muscle = "#34688d")

#------------------------
# Read in the objects
#------------------------

so.harmony <- readRDS(paste0(wdir, "20250708_10-26_Visualisation_.dir/RDS_objects.dir/Achilles_integrated_annotated.rds"))
so.achilles <- readRDS(paste0(wdir, "data/Achilles_fine_annotation.rds"))

# need to save and read in the subsetted objects.
# looks like I should also do the fine annotation on the subsets for muscle and stromal cells


#--------------------------------------------------
# Figure S2A: DimPlot of predicted broad annotations
#---------------------------------------------------

do_DimPlot(so.harmony, reduction = "harmony.umap", group.by = "predicted.id_integrated_broad",
           colors.use = achilles.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right"
)


ggsave(paste0(directory, "/Supplementary/Figure_S4/UMAP_predicted_broad.png"), width = 10, height = 7, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S4/UMAP_predicted_broad.svg"), width = 10, height = 7, bg = "white")


#--------------------------------------------------
# Figure S2B: Cell proportions across microanatomy
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

df <- left_join(df.sn, df.broad, by = c("broad_annotation" = "predicted.id_integrated_broad"))
df.number <- df %>% select(broad_annotation, snRNAseq, xenium)
df.number <- df.number %>% pivot_longer(2:3, names_to = "technology", values_to = "cell_number")

ggplot(df.number, aes(y=cell_number, x=broad_annotation, fill = technology))+
  geom_bar(position = "dodge", stat = "identity")+
  #geom_text(aes(label=cell_number),vjust=-1)+
  theme_classic()+
  scale_fill_manual(values = c("snRNAseq" = "#fc8d62", "xenium" = "#8da0cb"))+
  #ylim(c(0, 5500))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank())

ggsave((paste0(directory, "/Supplementary/Figure_S4/Cell_numbers.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S4/Cell_numbers.svg")), 
       width = 8, height = 5)

df.proportion <- df %>% select(broad_annotation, snRNAseq_percent, xenium_percent)
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

# Cell proportions for broad annotation split by microanatomy


df <- so.harmony[[]] %>% 
  select("predicted.id_integrated_broad", "microanatomical_site")%>% 
  group_by(predicted.id_integrated_broad, microanatomical_site) %>%  
  summarise(cell_number=n()) %>% 
  mutate(percent_cells=(cell_number/sum(cell_number))*100)


df$predicted.id_integrated_broad <- factor(df$predicted.id_integrated_broad, 
                              levels = df.broad$predicted.id_integrated_broad)

ggplot(df, aes(y=cell_number, x=predicted.id_integrated_broad, fill = microanatomical_site))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_manual(values = ma.cols)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank())

ggsave((paste0(directory, "/Supplementary/Figure_S4/Xenium_broad_annotation_microanatomy.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S4/Xenium_broad_annotation_microanatomy.svg")), 
       width = 8, height = 5)

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
# Figure S5D: Muscle fine annotation
#---------------------------------------------------

Idents(so.harmony) <- "predicted.id_integrated_fine"
so.musc <- subset(so.harmony, idents = c("Fast-twitch skeletal muscle cells", 
                                         "Slow-twitch skeletal muscle cells", 
                                         "Transitional skeletal muscle cells", 
                                         "Satellite cells"))

do_DimPlot(so.musc, reduction = "harmony.umap", group.by = "predicted.id_integrated_fine",
           #colors.use = muscle.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right"
)

ggsave(paste0(directory, "/Supplementary/Figure_S5/UMAP_predicted_muscle.png"), 
       width = 8, height = 4, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S5/UMAP_predicted_muscle.svg"), width = 10, height = 7, bg = "white")


#--------------------------------------------------
# Figure S6D: Stromal fine annotation
#---------------------------------------------------

# TODO

so.str <- subset(so.harmony, idents = c("Lymphatic endothelial cells", 
                                        "Adipocytes", 
                                        "vSMC", 
                                        "Arteriolar VEC",
                                        "Pericytes", 
                                        "Nervous system cells",
                                        "Venular VEC"))

do_DimPlot(so.str, reduction = "harmony.umap", group.by = "predicted.id_integrated_fine",
           #colors.use = muscle.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right"
)

ggsave(paste0(directory, "/Supplementary/Figure_S5/UMAP_predicted_stromal.png"), 
       width = 8, height = 4, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S5/UMAP_predicted_stromal.svg"), 
       width = 8, height = 4, bg = "white")



#--------------------------------------------------
# Figure S5D: Immune fine annotation
#---------------------------------------------------

# TODO

#--------------------------------------------------
# Figure S4D: Fibroblast fine annotation UMAP
#---------------------------------------------------

# UMAP to do

#--------------------------------------------------
# Figure S5G: Fibroblast cell proportions in snRNAseq and Xenium
#---------------------------------------------------

Idents(so.harmony) <- so.harmony$predicted.id_integrated_broad
so.fb <- subset(so.harmony, idents = "Fibroblasts")
df <- so.fb[[]] %>% select("predicted.id_integrated_broad", "fine_annotation")

df.fine <- df %>% 
  group_by(fine_annotation) %>%  
  summarise(xenium=n()) %>% 
  mutate(xenium_percent=(xenium/sum(xenium))*100) %>% 
  arrange(desc(xenium))

df.fine$fine_annotation <- factor(df.fine$fine_annotation, 
                                  levels = df.fine$fine_annotation)

df.sn <- so.achilles[[]] %>% select(fine_annotation)

df.sn <- df.sn %>% 
  group_by(fine_annotation) %>%  
  summarise(snRNAseq=n()) %>% 
  mutate(snRNAseq_percent=(snRNAseq/sum(snRNAseq))*100) %>% 
  arrange(desc(snRNAseq))
df.sn$fine_annotation <- factor(df.sn$fine_annotation, 
                                  levels = df.fine$fine_annotation)

df <- left_join(df.fine, df.sn)
df.number <- df %>% select(fine_annotation, snRNAseq, xenium)
df.number <- df.number %>% pivot_longer(2:3, names_to = "technology", values_to = "cell_number")

ggplot(df.number, aes(y=cell_number, x=fine_annotation, fill = technology))+
  geom_bar(position = "dodge", stat = "identity")+
  #geom_text(aes(label=cell_number),vjust=-1)+
  theme_classic()+
  scale_fill_manual(values = c("snRNAseq" = "#fc8d62", "xenium" = "#8da0cb"))+
  #ylim(c(0, 5500))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank())

ggsave((paste0(directory, "/Supplementary/Figure_S4/Fibroblast_cell_numbers.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S4/Fibroblast_cell_numbers.svg")), 
       width = 8, height = 5)

df.proportion <- df %>% select(fine_annotation, snRNAseq_percent, xenium_percent)
df.proportion <- df.proportion %>% pivot_longer(2:3, names_to = "technology", values_to = "percent_cells")

ggplot(df.proportion, aes(y=percent_cells, x=fine_annotation, fill = technology))+
  geom_bar(position = "dodge", stat = "identity")+
  theme_classic()+
  scale_fill_manual(values = c("snRNAseq_percent" = "#d95f02", "xenium_percent" = "#7570b3"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank())

ggsave((paste0(directory, "/Supplementary/Figure_S4/Fibroblast_cell_percent.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S4/Fibroblast_cell_percent.svg")), 
       width = 8, height = 5)

# Cell proportions for broad annotation split by microanatomy


df <- so.fb[[]] %>% 
  select("fine_annotation", "microanatomical_site")%>% 
  group_by(fine_annotation, microanatomical_site) %>%  
  summarise(cell_number=n()) %>% 
  mutate(percent_cells=(cell_number/sum(cell_number))*100)


df$fine_annotation <- factor(df$fine_annotation, 
                                           levels = df.fine$fine_annotation)

ggplot(df, aes(y=cell_number, x=fine_annotation, fill = microanatomical_site))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_manual(values = ma.cols)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank())

ggsave((paste0(directory, "/Supplementary/Figure_S4/Xenium_fibroblast_annotation_microanatomy.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S4/Xenium_fibroblast_annotation_microanatomy.svg")), 
       width = 8, height = 5)

ggplot(df, aes(y=percent_cells, x=fine_annotation, fill = microanatomical_site))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_fill_manual(values = ma.cols)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x = "", y = "Proportion")+
  labs(fill = "Microanatomical site")

ggsave((paste0(directory, "/Supplementary/Figure_S4/Xenium_fibroblast_microanatomy_percent.png")), 
       width = 8, height = 5)
ggsave((paste0(directory, "/Supplementary/Figure_S4/Xenium_fibroblast_microanatomy_percent.svg")), 
       width = 8, height = 5)

