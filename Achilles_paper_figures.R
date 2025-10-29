##############################
# Figures for Achilles paper #
##############################

#--------
# set up
#--------

library(Seurat)
library(tidyverse)
library(pals)
library("scales")
library("ggsci")
library(SCpubr)
library(miloR)
library(RColorBrewer)
library(scDotPlot)
library(SingleCellExperiment)
library(CrossTalkeR)
library(viridis)
library(liana)
library(data.table)

#------------------------
# make output directories
#------------------------

directory <- "Achilles_paper_figures.dir"
dir.create(directory)
dir.create(paste0(directory, "/Figure_1"))
dir.create(paste0(directory, "/Figure_1/UMAPs/"))
dir.create(paste0(directory, "/Figure_2"))
dir.create(paste0(directory, "/Figure_3"))
dir.create(paste0(directory, "/Figure_4"))
dir.create(paste0(directory, "/Supplementary"))
dir.create(paste0(directory, "/Supplementary/Heatmaps/"))
dir.create(paste0(directory, "/Supplementary/Annotation_dotplots/"))
dir.create(paste0(directory, "/Supplementary/Beeswarm_plots/"))
dir.create(paste0(directory, "/Supplementary/Figure_S1"))
dir.create(paste0(directory, "/Supplementary/Figure_S2/"))
dir.create(paste0(directory, "/Supplementary/Figure_S3/"))
dir.create(paste0(directory, "/Supplementary/Figure_S4/"))
dir.create(paste0(directory, "/Supplementary/Figure_S5/ma_specific_genes"))
dir.create(paste0(directory, "/Supplementary/Figure_S5/"))
dir.create(paste0(directory, "/Supplementary/Figure_S6/"))
dir.create(paste0(directory, "/Supplementary/Figure_S7/"))
dir.create(paste0(directory, "/Supplementary/Figure_S8/"))



#------------------------
# Read in the objects
#------------------------

# integrated with scVI, 5000 vf, decontX filter, soupX filter, fine annotation
so.achilles <- readRDS("/project/tendonhca/shared/chromium/analysis/20231006_achilles/20250228_12-50_Fine_annotation.dir/RDS_objects.dir/Achilles_fine_annotation.rds")

# fibroblast subset
so.fibroblast <- readRDS("20250228_12-50_Fine_annotation.dir/RDS_objects.dir/Achilles_fibroblast_subset.rds")

# muscle subset
so.muscle <- readRDS("20250228_12-50_Fine_annotation.dir/RDS_objects.dir/Achilles_muscle_subset.rds")

# immune subset
so.immune <- readRDS("20250228_12-50_Fine_annotation.dir/RDS_objects.dir/Achilles_immune_subset.rds")

# stromal subset
so.stromal <- readRDS("20250228_12-50_Fine_annotation.dir/RDS_objects.dir/Achilles_stromal_subset.rds")

# Simplify fibroblast cell names

df <- data.frame("fine_annotation_fibroblasts" = levels(so.fibroblast$fine_annotation_fibroblasts))
df$cell_annotation_pseudobulk <- c("NEGR1.VCAN.fb", "COMP.MMP3.fb", "NEGR1.COL15A1.fb", 
                                   "Chondrocytes", "NEGR1.ITGA6.fb", "COMP.THBS4.fb", "PRG4.fb")
metadata <- so.fibroblast[[]]
metadata <- metadata %>% left_join(df, by = "fine_annotation_fibroblasts")
so.fibroblast <- AddMetaData(so.fibroblast, metadata$cell_annotation_pseudobulk, "cell_annotation_pseudobulk")

#-------------------------------------------------
# Rename Nervous system cells as Schwann cells
#-------------------------------------------------

df <- so.achilles[[]] %>% select(broad_annotation, fine_annotation)
df$broad_annotation <- str_replace(df$broad_annotation, "Nervous system cells", "Schwann cells")
df$fine_annotation <- str_replace(df$fine_annotation, "Nervous system cells", "Schwann cells")
so.achilles$broad_annotation <- df$broad_annotation
so.achilles$fine_annotation <- df$fine_annotation

df <- so.stromal[[]] %>% select(fine_annotation_stromal)
df$fine_annotation_stromal <- str_replace(df$fine_annotation_stromal, "Nervous system cells", "Schwann cells")
so.stromal$fine_annotation_stromal <- df$fine_annotation_stromal
Idents(so.stromal) <- so.stromal$fine_annotation_stromal


#------------------------
# set the colours
#------------------------

### Using seaborn ###

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
                        "PRG4hi fibroblasts"= "#9267BF",
                        "Chondrocytes" = "#CD2321"
)

fibroblast.colours.short <- fibroblast.colours
names(fibroblast.colours.short) <- c("COMP.MMP3.fb",
                                     "COMP.THBS4.fb",
                                     "NEGR1.COL15A1.fb",
                                     "NEGR1.VCAN.fb", 
                                     "NEGR1.ITGA6.fb",
                                     "PRG4.fb",
                                     "Chondrocytes")



# Immune colours

immune.colours <- c("B cells" = "#88564A", 
                    "Monocytes" = "#21EFA5",
                    "CLEC9Ahi DCs" = "#3443eb", 
                    "MERKhi LYVE1lo macrophages" = "#C19C93", 
                    "T cells" = "#FA2CF6", 
                    "Slow twitch skeletal muscle cells" = "#F5F5F5",
                    "CLEC10Ahi DCs" = "#E8EF21",
                    "Granulocytes" = "#CD2321",
                    "MERTKhi LYVE1hi macrophages" = "#43BCE7",
                    "NK cells" = "#F9BA73",
                    "MERTKlo PTPRGhi macrophages" = "#003C30",
                    "Plasma cells" = "#9267BF"
)

# stromal colours

stromal.colours <- c("Lymphatic endothelial cells" = "#FE9A30", 
                     "Adipocytes" = "#D2C800",
                     "Arteriolar VEC" = "#F2B6D2",
                     "Venular VEC" =  "#DC77C3",
                     "vSMC" = "#9FDF86",
                     "Pericytes" = "#41A021",
                     "Schwann cells" = "#458BFF"
                     #"Nervous system cells" = "#F67E00"
                     )

# Muscle colours

muscle.colours <- c("Slow-twitch skeletal muscle cells" = "#A177FF",
                    "Transitional skeletal muscle cells" = "#39F0F0", 
                    "Satellite cells" = "#FF7B5A", 
                    "Fast-twitch skeletal muscle cells" = "#C2DD14")



# fine annotation colours
achilles.fine.colours <- c(fibroblast.colours, muscle.colours, stromal.colours, immune.colours)

# Microanatomy colours
#ma.cols <-  c(Enth = "#E98500", MB = "#009E73", MTJ = "#CC79A7", muscle = "#0072B2")
ma.cols <-  c(Enthesis = "#d37504", MB = "#959295", MTJ = "#359c72", Muscle = "#34688d")

# Xenium colours

xenium.explorer.colours <- c("#FFABC3",
                             "#FC6400", 
                             "#D2C800", 
                             "#7CEC37", 
                             "#21EFA5", 
                             "#43BCE7", 
                             "#458BFF", 
                             "#A177FF", 
                             "#DD51AD", 
                             "#ABABAB")

xenium.colours <- c("Fibroblasts" = "#FC6400", 
                    "COMPhi MMP3hi fibroblasts" = "#F9BA73", 
                    "NEGR1hi VCANhi fibroblasts" = "#FC6400",
                    "NEGR1hi ITGA6hi fibroblasts" = "#C19C93", 
                    "NEGR1hi COL15A1hi fibroblasts" = "#41A021",
                    "Skeletal muscle cells" = "#A177FF",
                    "Immune" = "#44BED0",
                    "Mural cells" = "#7CEC37",
                    "Adipocytes" = "#D2C800",
                    "Endothelial cells" = "#FFABC3",
                    "Schwann cells" = "#458BFF", 
                    "Satellite cells" = "#C4B0D6",
                    "Macrophages" = "#44BED0", 
                    "B cells" = "#88564A", 
                    "T cells" = "#BBBD00", 
                    "Plasma cells" = "#01EBD6", 
                    "Granulocytes" = "#C7C7C7", 
                    # Stromal
                    "Vascular endothelial cells" = "#FFABC3",
                    "Lymphatic endothelial cells" = "#DD51AD"
                    )



#####################################
# Figure 1 - overview of cell types #
#####################################

#-----------------------------------------------------------------------
# Figure 1A: # Diagram of Achilles tendon microanatomy 
#-----------------------------------------------------------------------

#----------------------------------------------
# Figure 1B: UMAP of Achilles broad annotation
#----------------------------------------------

Idents(so.achilles) <- so.achilles$broad_annotation
do_DimPlot(so.achilles, reduction = "umap.scvi", group.by = "broad_annotation",
           colors.use = achilles.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "none"
           )
ggsave(paste0(directory, "/Figure_1/UMAP_broad.png"), width = 10, height = 8, bg = "white")
ggsave(paste0(directory, "/Figure_1/UMAP_broad.svg"), width = 10, height = 8, bg = "white")

#--------------------------------------------
# Figure 1C: MiloR results (broad annotation)
#--------------------------------------------

da_results <- read.table("20250226_10-04_MiloR_Cell_Proportions.dir/Results.dir/Achilles-milo-DA-results.txt", 
                         sep = "\t", header = TRUE)
da_results <- da_results %>% filter(broad_annotation != "Mixed")
da_results$broad_annotation <- str_replace(da_results$broad_annotation, "Nervous system cells", "Schwann cells")

da_results$broad_annotation <- factor(da_results$broad_annotation, 
                                      levels = c("Schwann cells", 
                                                 "Adipocytes", 
                                                 "Mural cells", 
                                                 "Lymphatic endothelial cells", 
                                                 "Vascular endothelial cells", 
                                                 "Plasma cells",
                                                 "B cells", 
                                                 "T cells", 
                                                 "Granulocytes", 
                                                 "Macrophages", 
                                                 "Satellite cells", 
                                                 "Skeletal muscle cells", 
                                                 "Fibroblasts"))
plotDAbeeswarm(da_results, group.by = "broad_annotation")+
    theme(axis.title.y=element_blank())+
    scale_fill_gradient2(low = "#00466B", mid = "#FFFFFF", high = "#A91400", midpoint = 0, name = "logFC")
ggsave(paste0(directory, "/Figure_1/Beeswarm_plot_broad.png"), width = 10, height = 10, bg = "white")
ggsave(paste0(directory, "/Figure_1/Beeswarm_plot_broad.svg"), width = 10, height = 10, bg = "white")

#-------------------------------------------------------------------
# Figure 1D: UMAP of Achilles broad annotation split by microanatomy
#-------------------------------------------------------------------


# this plots all the microanatomies but it is too big to deal with in Affinity Designer as svg
# use .png instead
do_DimPlot(so.achilles, reduction = "umap.scvi", group.by = "broad_annotation", 
           split.by = "microanatomical_site", 
           colors.use = achilles.colours, 
           plot_cell_borders = TRUE, 
           pt.size = 0.5, 
           legend.position = "right", 
           ncol = 5, 
           ) # font.type "serif" is times new roman, "mono" is courier new
ggsave(paste0(directory, "/Figure_1/UMAPs/UMAP_by_microanatomy.png"), width = 32, height = 6, bg = "white")
#ggsave(paste0(directory, "/Figure_1/UMAPs/UMAP_by_microanatomy.svg"), width = 25, height = 10, bg = "white")


do_DimPlot(so.achilles, reduction = "umap.scvi", group.by = "broad_annotation", 
           split.by = "microanatomical_site", 
           colors.use = xenium.colours, 
           plot_cell_borders = TRUE, 
           pt.size = 0.5, 
           legend.position = "none", 
           ncol = 1, 
) # font.type "serif" is times new roman, "mono" is courier new
ggsave(paste0(directory, "/Figure_1/UMAPs/UMAP_by_microanatomy_tall.png"), width = 8, height = 25, bg = "white")
#ggsave(paste0(directory, "/Figure_1/UMAPs/UMAP_by_microanatomy_tall.svg"), width = 6, height = 25, bg = "white")

# NB cannot get rid of "original" unless use new github version of package
# https://github.com/enblacar/SCpubr/issues/45
# Also the spacing of the titles is wrong

#-------------------------------------------------------------------
# Figure 1E, F: Xenium cell types & H&E, made in Xenium Explorer
#-------------------------------------------------------------------

#############################
# Figure 2: Fibroblasts #
#############################

#----------------------------------------------------
# Figure S2A: Clustered gene expression in fibroblasts
#----------------------------------------------------

cluster.fb <- readRDS(paste0("20250304_17-29_Pseudobulk.dir/Cluster_results/Fibroblasts_cluster_results.rds"))
cluster.fb$plot+
    facet_wrap("cluster", ncol = 5)+
    theme_classic()+
    theme(legend.position = "none")+
    scale_colour_manual(values = "#00466B")
ggsave(paste0(directory, "/Figure_2/Fibroblasts_clusters.png"), width = 16, height = 8)
ggsave(paste0(directory, "/Figure_2/Fibroblasts_clusters.svg"), width = 16, height = 8)

#------------------------------------------------------------------------
# Figure S2B: Pathway analysis on clustered gene expression in fibroblasts
#------------------------------------------------------------------------

# read in combined results from Fibroblasts

pathways.fb <- read.table("20250305_10-46_Pathway_Analysis.dir/Results/Fibroblast_results/Fibroblasts_up_vs_down.txt", 
                          sep = "\t", header = TRUE)

# make comparison  a factor
pathways.fb$comparison <- factor(pathways.fb$comparison, levels = c("upregulated", "downregulated", "tendon_high"))

# make IDs a factor
pathways <- as.data.frame(sort(table(pathways.fb$ID)))
colnames(pathways) <- c("ID", "pathway_freq")
pathways.fb <- pathways.fb %>% 
    left_join(pathways) %>%
    group_by(pathway_freq) %>%
    arrange(comparison)
# here could add a filter for the min number of genes overlapping with the pathway

pathways.fb$ID <- factor(pathways.fb$ID, levels = rev(unique(pathways.fb$ID)))

ggplot(pathways.fb, aes(x = comparison, y = ID, size = Count, colour = new_col)) +
    geom_point()+
    #scale_colour_viridis()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_blank())+
    labs(colour = "-log10(p.adj)")+
    scale_y_discrete(position = "right")

ggsave (paste0(directory, "/Figure_2/Fibroblasts_up_vs_down.png"), 
        width =6, height = 7)
ggsave (paste0(directory, "/Figure_2/Fibroblasts_up_vs_down.svg"), 
        width = 6, height = 7)



#---------------------------------------------
# Figure 2C: Fibroblasts split by microanatomy
#---------------------------------------------

do_DimPlot(so.fibroblast, reduction = "umap", group.by = "fine_annotation_fibroblasts",
           split.by = "microanatomical_site",
           colors.use = fibroblast.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5,
           ncol = 5,
           legend.position = "right")
ggsave(paste0(directory, "/Figure_2/UMAP_fibroblasts_by_microanatomy.png"), width = 25, height = 5,  bg = "white")
ggsave(paste0(directory, "/Figure_2/UMAP_fibroblasts_by_microanatomy.svg"), width = 25, height = 5, , bg = "white")

#----------------------------------------------
# Figure 2D: MiloR beeswarm plot on fibroblasts
#----------------------------------------------

# MiloR on fibroblast subset

da_results <- read.table("20250226_11-09_MiloR_Cell_Proportions.dir/Results.dir/Achilles-milo-DA-results.txt", 
                         sep = "\t", header = TRUE)
da_results <- da_results %>% filter(fine_annotation_fibroblasts != "Mixed")
da_results$fine_annotation_fibroblasts <- factor(da_results$fine_annotation_fibroblasts,
                                                 levels = c("NEGR1hi COL15A1hi fibroblasts", 
                                                            "COMPhi THBS4hi fibroblasts",     
                                                            "NEGR1hi ITGA6hi fibroblasts",   
                                                            "PRG4hi fibroblasts",
                                                            "NEGR1hi VCANhi fibroblasts",
                                                            "COMPhi MMP3hi fibroblasts",
                                                            "Chondrocytes"))

plotDAbeeswarm(da_results, group.by = "fine_annotation_fibroblasts")+
    theme(axis.title.y=element_blank())+
    scale_fill_gradient2(low = "#00466B", mid = "#FFFFFF", high = "#A91400", midpoint = 0, name = "logFC")
ggsave(paste0(directory, "/Figure_2/Beeswarm_plot_fibroblasts.png"), width = 10, height = 8, bg = "white")
ggsave(paste0(directory, "/Figure_2/Beeswarm_plot_fibroblasts.svg"), width = 10, height = 8, bg = "white")






#---------------------------------------------------
# Figure 2E: Xenium and H&E made in Xenium explorer
#---------------------------------------------------

#---------------------------------------------------
# Figure 2F: COMP and NEGR1 across spatial
#---------------------------------------------------

# Xenium script

#############################################
# Figure 3: Fine annotation #
#############################################

#------------------------------------
# Figure 3A: Stromal UMAP & Beeswarm
#------------------------------------

do_DimPlot(so.stromal, reduction = "umap", group.by = "fine_annotation_stromal",
           colors.use = stromal.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right")
ggsave(paste0(directory, "/Figure_3/UMAP_stromal.png"), width = 10, height = 7, bg = "white")
ggsave(paste0(directory, "/Figure_3/UMAP_stromal.svg"), width = 10, height = 7, bg = "white")


# MiloR on stromal subset

da_results <- read.table("20250226_12-47_MiloR_Cell_Proportions.dir/Results.dir/Achilles-milo-DA-results.txt", 
                         sep = "\t", header = TRUE)
da_results <- da_results %>% filter(fine_annotation_stromal != "Mixed")
da_results$fine_annotation_stromal <- factor(da_results$fine_annotation_stromal,
                                             levels = c("Pericytes", 
                                                        "Lymphatic endothelial cells", 
                                                        "Schwann cells",
                                                        "Adipocytes", 
                                                        "vSMC", 
                                                        "Arteriolar VEC", 
                                                        "Venular VEC"
                                             ))

plotDAbeeswarm(da_results, group.by = "fine_annotation_stromal")+
    theme(axis.title.y=element_blank())+
    scale_fill_gradient2(low = "#00466B", mid = "#FFFFFF", high = "#A91400", midpoint = 0, name = "logFC")+
    ylim(c(-4, 4))
ggsave(paste0(directory, "/Figure_3/Beeswarm_plots/Beeswarm_plot_stromal.png"), width = 10, height = 8, bg = "white")
#ggsave(paste0(directory, "/Figure_3/Beeswarm_plots/Beeswarm_plot_stromal.svg"), width = 10, height = 8, bg = "white")


#-----------------------------------
# Figure 3B: Immune UMAP & Beeswarm
#-----------------------------------

# remove muscle cells for plotting

so.immune.subset <- subset(so.immune, idents = "Slow-twitch skeletal muscle cells", invert = TRUE)

do_DimPlot(so.immune.subset, reduction = "umap", group.by = "fine_annotation_immune",
           colors.use = immune.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right")
ggsave(paste0(directory, "/Figure_3/UMAP_immune.png"), width = 10, height = 7, bg = "white")
ggsave(paste0(directory, "/Figure_3/UMAP_immune.svg"), width = 10, height = 7, bg = "white")

# MiloR on immune subset
da_results <- read.table("20250226_12-02_MiloR_Cell_Proportions.dir/Results.dir/Achilles-milo-DA-results.txt", 
                         sep = "\t", header = TRUE)
da_results <- da_results %>% filter(fine_annotation_immune != "Mixed")
da_results$fine_annotation_immune <- factor(da_results$fine_annotation_immune,
                                            levels = c("MERTKhi LYVE1hi macrophages", 
                                                       "T cells",
                                                       "Granulocytes", 
                                                       "MERKhi LYVE1lo macrophages", 
                                                       "Monocytes", 
                                                       "Plasma cells",
                                                       "CLEC10A DCs",
                                                       "B cells",
                                                       "MERTKlo PTPRGhi macrophages"))

plotDAbeeswarm(da_results, group.by = "fine_annotation_immune")+
    theme(axis.title.y=element_blank())+
    scale_fill_gradient2(low = "#00466B", mid = "#FFFFFF", high = "#A91400", midpoint = 0, name = "logFC")+
    ylim(c(-3, 3))
dir.create(paste0(directory, "/Figure_3/Beeswarm_plots/"))
ggsave(paste0(directory, "/Figure_3/Beeswarm_plots/Beeswarm_plot_immune.png"), width = 10, height = 8, bg = "white")
#ggsave(paste0(directory, "/Figure_3/Beeswarm_plots/Beeswarm_plot_immune.svg"), width = 10, height = 8, bg = "white")

#----------------------------------
# Figure 3C: Muscle UMAP & Beeswarm
#----------------------------------

do_DimPlot(so.muscle, reduction = "umap", group.by = "fine_annotation_muscle",
           colors.use = muscle.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right")
ggsave(paste0(directory, "/Figure_3/UMAP_muscle.png"), width = 10, height = 7, bg = "white")
ggsave(paste0(directory, "/Figure_3/UMAP_muscle.svg"), width = 10, height = 7, bg = "white")

do_DimPlot(so.muscle, reduction = "umap", group.by = "fine_annotation_muscle",
           split.by = "microanatomical_site",
           colors.use = muscle.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right")

ggsave(paste0(directory, "/Figure_3/UMAP_muscle_microanatomy.png"), width = 12, height = 7, bg = "white")

# MiloR on muscle subset
# NB this does not include the enthesis samples 

da_results <- read.table("20250226_15-10_MiloR_Cell_Proportions.dir/Results.dir/Achilles-milo-DA-results.txt", 
                         sep = "\t", header = TRUE)
da_results <- da_results %>% filter(fine_annotation_muscle != "Mixed")
da_results$fine_annotation_muscle <- factor(da_results$fine_annotation_muscle,
                                            levels = c("Fast-twitch",
                                                       "Satellite cells",
                                                       "Slow-twitch",
                                                       "MTJ-specific"))

plotDAbeeswarm(da_results, group.by = "fine_annotation_muscle")+
    theme(axis.title.y=element_blank())+
    scale_fill_gradient2(low = "#00466B", mid = "#FFFFFF", high = "#A91400", midpoint = 0, name = "logFC")+
    ylim(c(-6, 6.5))
ggsave(paste0(directory, "/Figure_3/Beeswarm_plots/Beeswarm_plot_muscle.png"), width = 10, height = 8, bg = "white")
#ggsave(paste0(directory, "/Figure_3/Beeswarm_plots/Beeswarm_plot_muscle.svg"), width = 10, height = , bg = "white")

#-----------------------------
# Figure 3D-F: Xenium Explorer
#-----------------------------



#####################################
# Figure 4: Cell-cell communication #
#####################################

#--------------------------------------------------------
# Figure 4A: Cell-cell communication across microanatomy
#--------------------------------------------------------

# update some colours that are a bit similar
achilles.fine.colours["Pericytes"] <- "#2c6129"
achilles.fine.colours["MERTKhi LYVE1hi macrophages"] <- "#222b82"
achilles.fine.colours["MERKhi LYVE1lo macrophages"] <- "#faf89b"
achilles.fine.colours["NK cells"] <- "#964d03"
achilles.fine.colours["Slow twitch skeletal muscle cells"] <- "#A177FF"

crosstalker_out <- readRDS("20250311_14-50_Liana-fine_annotation.dir/Results/Crosstalker_report.rds")

# Enthesis
png(file=paste0(directory, "/Figure_4/Enthesis.png"), width=1600, height=1000)
p1 <- plot_cci(graph = crosstalker_out@graphs$Enth,
               colors = achilles.fine.colours, # can provide a named vector of custom colours
               plt_name = 'Enthesis',
               coords = crosstalker_out@coords[igraph::V(crosstalker_out@graphs$Enth)$name,],
               emax = 65,
               leg = TRUE,
               log = FALSE,
               pg = crosstalker_out@rankings[["Enth"]]$Pagerank,
               vnames = FALSE)
dev.off()

svg(file=paste0(directory, "/Figure_4/Enthesis.svg"), width=1600, height=1000)
p1 <- plot_cci(graph = crosstalker_out@graphs$Enth,
         colors = achilles.fine.colours, # can provide a named vector of custom colours
         plt_name = 'Enthesis',
         coords = crosstalker_out@coords[igraph::V(crosstalker_out@graphs$Enth)$name,],
         emax = 65,
         leg = TRUE,
         log = FALSE,
         pg = crosstalker_out@rankings[["Enth"]]$Pagerank,
         vnames = FALSE)
dev.off()

# Midbody

png(file=paste0(directory, "/Figure_4/Midbody.png"), width=1600, height=1000)
p2 <- plot_cci(graph = crosstalker_out@graphs$MB,
               colors = achilles.fine.colours,
               plt_name = 'Midbody',
               coords = crosstalker_out@coords[igraph::V(crosstalker_out@graphs$MB)$name,],
               emax = 65,
               leg = TRUE,
               pg = crosstalker_out@rankings[["MB"]]$Pagerank,
               vnames = FALSE)

dev.off()

svg(file=paste0(directory, "/Figure_4/Midbody.svg"), width=1600, height=1000)
p2 <- plot_cci(graph = crosstalker_out@graphs$MB,
               colors = achilles.fine.colours,
               plt_name = 'Midbody',
               coords = crosstalker_out@coords[igraph::V(crosstalker_out@graphs$MB)$name,],
               emax = 65,
               leg = TRUE,
               pg = crosstalker_out@rankings[["MB"]]$Pagerank,
               vnames = FALSE)

dev.off()

# MTJ

png(file=paste0(directory, "/Figure_4/MTJ.png"), width=1600, height=1000)
p3 <- plot_cci(graph = crosstalker_out@graphs$MTJ,
               colors = achilles.fine.colours,
               plt_name = 'MTJ',
               coords = crosstalker_out@coords[igraph::V(crosstalker_out@graphs$MTJ)$name,],
               emax = 65,
               leg = TRUE,
               pg = crosstalker_out@rankings[["MTJ"]]$Pagerank,
               vnames = FALSE)
dev.off()

svg(file=paste0(directory, "/Figure_4/MTJ.svg"), width=1600, height=1000)
p3 <- plot_cci(graph = crosstalker_out@graphs$MTJ,
               colors = achilles.fine.colours,
               plt_name = 'MTJ',
               coords = crosstalker_out@coords[igraph::V(crosstalker_out@graphs$MTJ)$name,],
               emax = 65,
               leg = TRUE,
               pg = crosstalker_out@rankings[["MTJ"]]$Pagerank,
               vnames = FALSE)
dev.off()

# Muscle

png(file=paste0(directory, "/Figure_4/Muscle.png"), width=1600, height=1000)
p4 <- plot_cci(graph = crosstalker_out@graphs$muscle,
               colors = achilles.fine.colours,
               plt_name = 'Muscle',
               coords = crosstalker_out@coords[igraph::V(crosstalker_out@graphs$muscle)$name,],
               emax = 65,
               leg = TRUE,
               pg = crosstalker_out@rankings[["muscle"]]$Pagerank,
               vnames = TRUE)
dev.off()

svg(file=paste0(directory, "/Figure_4/Muscle.svg"), width=1600, height=1000)
p4 <- plot_cci(graph = crosstalker_out@graphs$muscle,
               colors = achilles.fine.colours,
               plt_name = 'Muscle',
               coords = crosstalker_out@coords[igraph::V(crosstalker_out@graphs$muscle)$name,],
               emax = 65,
               leg = TRUE,
               pg = crosstalker_out@rankings[["muscle"]]$Pagerank, 
               vnames = FALSE)
dev.off()

#-----------------------------
# Figure 4B: VEC interactions
#-----------------------------

liana_res <- read.table("20250903_10-53_Liana.dir/Liana_aggregate_results_fibroblast_annotation.txt", header = TRUE, sep = "\t")

fibroblasts <- c("COMP.MMP3.fb", 
                 "COMP.THBS4.fb", "NEGR1.VCAN.fb", 
                 "NEGR1.COL15A1.fb", 
                 "NEGR1.ITGA6.fb", 
                 "PRG4.fb", "Chondrocytes")
stromal <- c("VEC", "Mural cells", "Adipocytes", "LEC", "Schwann cells")
immune <- c("Macrophages", "NK cells", "Granulocytes", "B cells", "T cells")
muscle <- c("Skeletal muscle cells", "Satellite cells")
fb.colours <- c(fibroblast.colours.short[2:4], fibroblast.colours.short[6])
liana_res %>%
    filter(natmi.edge_specificity>0.1) %>% # add a specificity filter
    liana_dotplot(source_groups = "VEC",
                  target_groups = fibroblasts,
                  ntop = 50)+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, 
                                     colour = fb.colours))

ggsave(paste0(directory, "/Figure_4/VEC-fibroblast.png"), device = "png", 
       bg = "white", width = 11, height = 12)

ggsave(paste0(directory, "/Figure_4/VEC-fibroblast.svg"), device = "svg", 
       bg = "white", width = 12, height = 10)

#----------------------------------------
# Figure 4C: NEGR1 COL15A1 fb interactions
#----------------------------------------

str.cols <- c(stromal.colours[1], stromal.colours[7], stromal.colours[4])
liana_res %>%
    filter(natmi.edge_specificity>0.1) %>% # add a specificity filter
    liana_dotplot(source_groups = "NEGR1.COL15A1.fb",
                  target_groups = stromal,
                  ntop = 50)+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, 
                                     colour = str.cols))

ggsave(paste0(directory, "/Figure_4/COL15A1fb-stromal.png"), device = "png", 
       bg = "white", width = 11, height = 12)

ggsave(paste0(directory, "/Figure_4/COL15A1fb-stromal"), device = "svg", 
       bg = "white", width = 12, height = 10)


#----------------------------------------
# Figure 4D: Schwann cells interactions
#----------------------------------------
fb.colours <- c(fibroblast.colours.short[7], fibroblast.colours.short[1], 
                fibroblast.colours.short[3], 
                fibroblast.colours.short[5:6])

liana_res %>%
    filter(natmi.edge_specificity>0.1) %>% # add a specificity filter
    liana_dotplot(source_groups = "Schwann cells",
                  target_groups = fibroblasts,
                  ntop = 50)+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, colour = fb.colours))


ggsave(paste0(directory, "/Figure_4/Schwann-fibroblast.png"), device = "png", 
       bg = "white", width = 11, height = 12)

ggsave(paste0(directory, "/Figure_4/Schwann-fibroblast"), device = "svg", 
       bg = "white", width = 12, height = 10)

##########################
# Supplementary Figures  #
##########################


#----------------------
# Figure S1: H&E images
#----------------------

#-----------------------------------------
# Figure S2A: Quality control of snRNA-seq
#-----------------------------------------

do_FeaturePlot(so.achilles, reduction = "umap.scvi", features = c("decontX_contamination", "soupX_fraction"),
               plot_cell_borders = TRUE, pt.size = 0.5)
ggsave(paste0(directory, "/Supplementary/Figure_S2/UMAP_decontX_soupX.png"), width = 12, height = 7)


#---------------------------------
# Figure S2B: Integration by patient 
#---------------------------------

do_DimPlot(so.achilles, reduction = "umap.scvi", group.by = "broad_annotation",
           split.by = "patient",
           colors.use = achilles.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "none")
ggsave(paste0(directory, "/Supplementary/Figure_S2/UMAP_by_patient.png"), width = 12, height = 8)


#---------------------------
# Figure S2 and S6-8: Heatmaps
#---------------------------

# read in lists of top markers per cell type

broad_markers <- read.table("20250228_12-50_Fine_annotation.dir/Marker_lists/broad_annotation.txt", 
                            header = TRUE, sep = "\t")
fine_markers <- read.table("20250228_12-50_Fine_annotation.dir/Marker_lists/fine_annotation.txt", 
                           header = TRUE, sep = "\t")
fibroblast_markers <- read.table("20250228_12-50_Fine_annotation.dir/Marker_lists/fine_annotation_fibroblasts.txt", 
                                 header = TRUE, sep = "\t")
immune_markers <- read.table("20250228_12-50_Fine_annotation.dir/Marker_lists/fine_annotation_immune.txt", 
                             header = TRUE, sep = "\t")
muscle_markers <- read.table("20250228_12-50_Fine_annotation.dir/Marker_lists/fine_annotation_muscle.txt", 
                             header = TRUE, sep = "\t")
stromal_markers <- read.table("20250228_12-50_Fine_annotation.dir/Marker_lists/fine_annotation_stromal.txt", 
                              header = TRUE, sep = "\t")
# remove muscle markers from immune
immune_markers <- immune_markers %>% filter(cluster != "Slow-twitch skeletal muscle cells")

# change to Schwann cells
broad_markers$cluster <- str_replace(broad_markers$cluster, "Nervous system cells", "Schwann cells")
fine_markers$cluster <- str_replace(fine_markers$cluster, "Nervous system cells", "Schwann cells")
stromal_markers$cluster <- str_replace(stromal_markers$cluster, "Nervous system cells", "Schwann cells")


# plot heatmaps

do_heatmap <- function(so, marker_list, name, fig_width = 8, fig_height = 7){
    
    # make a factor of the cell names in order of descending frequency
    df <- as.data.frame(table(Idents(so)))%>% arrange(desc(Freq))
    df$Var1 <- factor(df$Var1, levels = df$Var1)
    
    # set the factor levels of the so
    levels(so) <- levels(df$Var1)
    
    # make the order of clusters a factor to match the Idents
    marker_list$cluster <- factor(marker_list$cluster, 
                                  levels = levels(Idents(so)))
    # select the colours
    colour_name <- get(paste0(name, ".colours"))
    print(colour_name)
    
    # filter top 10 genes
    marker_list_top10 <- marker_list %>%
        group_by(cluster) %>% 
        slice_max(n=10, order_by = abs(avg_log2FC)) %>% 
        pull(gene)
    
    # plot heatmap
    DoHeatmap(so, 
              features = unique(marker_list_top10),
              assay = "soupX", 
              group.colors = colour_name, 
              label = FALSE) + 
        scale_fill_viridis()
    ggsave(paste0(directory, "/Supplementary/Heatmaps/", name, "_heatmap.png"), width = fig_width, height = fig_height)
    #ggsave(paste0(directory, "/Supplementary/Heatmaps/", name, "_heatmap.svg"), width = fig_width, height = fig_width)
}

Idents(so.achilles) <- so.achilles$broad_annotation
do_heatmap(so.achilles, broad_markers, "achilles", fig_width = 12, fig_height = 15)
do_heatmap(so.fibroblast, fibroblast_markers, "fibroblast")
do_heatmap(so.immune.subset, immune_markers, "immune", fig_width = 12, fig_height = 12)
do_heatmap(so.muscle, muscle_markers, "muscle")
do_heatmap(so.stromal, stromal_markers, "stromal")


#----------------------------------------
# Figure S2B and S6-8: Clustered dotplots
#----------------------------------------

so_list <- list(so.stromal, so.immune, so.fibroblast, so.muscle, so.achilles)
sce_list <- lapply(so_list, as.SingleCellExperiment)
names(sce_list) <- c("stromal", "immune","fibroblast", "muscle", "achilles")
markers_list <- list(stromal_markers, immune_markers, fibroblast_markers, muscle_markers, broad_markers)
names(markers_list) <- c("stromal", "immune","fibroblast", "muscle", "achilles")

select_markers <- function(marker_list){
    # select the top 5 markers for each cluster
    features_5 <- marker_list %>% 
        group_by(cluster) %>%
        dplyr::slice(1:5) %>%
        dplyr::select(cluster, gene)
    
    # remove duplicated rows
    features_5 <- features_5[!duplicated(features_5$gene),]
    
    # make features into a named vector
    features_5 <- features_5 %>%
        deframe()
}

markers_list_top5 <- lapply(markers_list, select_markers)
names(markers_list_top5) <- c("stromal", "immune","fibroblast", "muscle", "achilles")


for (i in seq(1:length(sce_list))){
    rowData(sce_list[[i]])$Annotation <- markers_list_top5[[i]][match(rownames(sce_list[[i]]), markers_list_top5[[i]])] %>% 
        names()
}

# plot the top 5 markers (scaled data)

# broad annotation
sce_list[["achilles"]] %>%
    scDotPlot(features = markers_list_top5[["achilles"]],
              scale = TRUE,
              group = "broad_annotation",
              groupAnno = "broad_annotation",
              featureAnno = "Annotation",
              groupLegends = FALSE,
              annoColors = list("broad_annotation" = achilles.colours,
                                "Annotation" = achilles.colours),
              annoWidth = 0.02)
ggsave(paste0(directory, "/Supplementary/Figure_S2/broad_annotation_top5_dotplot.png"), 
       width = 10, height = 15, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S2/broad_annotation_top5_dotplot.svg"), 
       width = 10, height = 15, bg = "white")

# fibroblast
sce_list[["fibroblast"]] %>%
    scDotPlot(features = markers_list_top5[["fibroblast"]],
              scale = TRUE,
              group = "fine_annotation_fibroblasts",
              groupAnno = "fine_annotation_fibroblasts",
              featureAnno = "Annotation",
              groupLegends = FALSE,
              annoColors = list("fine_annotation_fibroblasts" = fibroblast.colours,
                                "Annotation" = fibroblast.colours),
              annoWidth = 0.02)
ggsave(paste0(directory, "/Supplementary/Figure_S5/fibroblast_top5_dotplot.png"), 
       width = 10, height = 10, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S5/fibroblast_top5_dotplot.svg"), 
       width = 10, height = 10, bg = "white")

# stromal
sce_list[["stromal"]] %>%
    scDotPlot(features = markers_list_top5[["stromal"]],
              scale = TRUE,
              group = "fine_annotation_stromal",
              groupAnno = "fine_annotation_stromal",
              featureAnno = "Annotation",
              groupLegends = FALSE,
              annoColors = list("fine_annotation_stromal" = stromal.colours,
                                "Annotation" = stromal.colours),
              annoWidth = 0.02)
ggsave(paste0(directory, "/Supplementary/Figure_S6/stromal_top5_dotplot.png"), 
       width = 10, height = 10, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S6/stromal_top5_dotplot.svg"), 
       width = 10, height = 10, bg = "white")

# immune
sce_list[["immune"]] %>%
    scDotPlot(features = markers_list_top5[["immune"]],
              scale = TRUE,
              group = "fine_annotation_immune",
              groupAnno = "fine_annotation_immune",
              featureAnno = "Annotation",
              groupLegends = FALSE,
              annoColors = list("fine_annotation_immune" = immune.colours,
                                "Annotation" = immune.colours),
              annoWidth = 0.02)
ggsave(paste0(directory, "/Supplementary/Figure_S7/immune_top5_dotplot.png"), 
       width = 10, height = 10, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S7/immune_top5_dotplot.svg"), 
       width = 10, height = 10, bg = "white")


# muscle
sce_list[["muscle"]] %>%
    scDotPlot(features = markers_list_top5[["muscle"]],
              scale = TRUE,
              group = "fine_annotation_muscle",
              groupAnno = "fine_annotation_muscle",
              featureAnno = "Annotation",
              groupLegends = FALSE,
              annoColors = list("fine_annotation_muscle" = muscle.colours,
                                "Annotation" = muscle.colours),
              annoWidth = 0.02)
ggsave(paste0(directory, "/Supplementary/Figure_S8/muscle_top5_dotplot.png"), 
       width = 10, height = 10, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S8/muscle_top5_dotplot.svg"), 
       width = 10, height = 10, bg = "white")



#---------------------------------------
# Figure S3C, S5C and S6-8: Annotation dotplots
#---------------------------------------

annotation_dotplot <- function (so, genes, category, fig_width = 7, fig_height = 7){
    markers <- unlist(genes, use.names = FALSE)
    p <- DotPlot(so, features = markers, assay = "soupX") + 
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "#f7f3f2", mid = "#f2390a", high = "#6e1e0a", midpoint = 1)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16))+
        theme(axis.text.y = element_text(size = 16))+
        theme(axis.title.x = element_blank())+
        theme(axis.title.y = element_blank())+
        theme(strip.text = element_blank())
    
    ggsave(paste0(directory, "/Supplementary/Annotation_dotplots/", category,".png"), bg = "white", 
           width = fig_width, height = fig_height)
    ggsave(paste0(directory, "/Supplementary/Annotation_dotplots/", category,".svg"), bg = "white", 
           width = fig_width, height = fig_height)
}

general_markers_top5 <- list("Fibroblast" = c("COL1A2", "COL3A1", "DCN", "NEGR1", "COMP"),
                             "Sk_musc" = c("TRDN", "TTN", "NEB", "DES", "ACTN2"), 
                             "VEC" = c("PECAM1", "PTPRB", "FLT1", "VWF", "ADGRL4"), 
                             "Macrophage"  = c("F13A1", "CD163", "MRC1", "MSR1","MERTK"), 
                             "Mural" = c( "NOTCH3", "PDGFRB", "MYO1B", "ELN", "FBN1"), 
                             "Adipocyte"  = c ("PLIN1", "LPL", "GPAM", "AQP7","ADIPOQ"),
                             "T cell"  = c( "CD247", "SKAP1", "THEMIS", "PTPRC", "IL7R"), 
                             "LEC" = c( "MMRN1", "PROX1", "PKHD1L1", "FLT4", "NRG3"),
                             "Satellite cells" = c("PAX7", "CALCR", "CDH4", "CLCN5", "CTNND2"),
                             "Granulocyte" = c("KIT", "CPA3", "IL18R1","CDK15", "MS4A2"), 
                             "B cell" = c("MS4A1", "CD37", "BLNK", "FCRL1", "IGHM"), 
                             "Schwann" = c("NRXN1", "XKR4","BMS1P14","CADM2","IL1RAPL2"),
                             "Plasma" = c("IGHG1", "MZB1", "XBP1", "PRDM1", "JCHAIN")
)


annotation_dotplot(so.achilles, general_markers_top5, "General_markers", 20, 8)


immune_markers <- list("T cell" = c("CD247", "SKAP1", "THEMIS", "PTPRC", "IL7R"), 
                       "MERTK+LYVE1-mac" = c("CTSL", "CTSB", "TPRG1", "HMOX1","ELL2"),
                       "MERTK-PTPRG+mac" = c("PTPRG", "MIR99AHG", "PARD3",
                                             "KAZN",  "IFI44L"), 
                       "MERTK+LYVE1+mac" = c("CD163", "F13A1","STAB1",
                                             "MERTK", "LYVE1"), 
                       "CLEC10A DC" = c("CLEC10A", "CD1C", "FCER1A", "IL1R2"),
                       "NK cell" = c("NCAM1", "GNLY", "KLRD1", "CCL5", "MCTP2"),
                       "Monocytes" = c("FCN1", "VCAN", "FCGR3A", "LST1", "LILRB2"),
                       "Granulocyte" = c( "KIT", "CPA3", "IL18R1","CDK15", "MS4A2"),
                       "Plasma" = c("IGHG1", "MZB1", "XBP1", "PRDM1", "JCHAIN"),
                       "CLEC9Ahi DC" = c("CLEC9A", "CADM1", "IDO1", "CLNK", "ZNF366"),
                       "B cells" = c("MS4A1", "CD37", "BLNK", "FCRL1", "IGHM")
                       
)
annotation_dotplot(so.immune.subset, immune_markers, "Immune_markers", 22, 8)

stromal_markers_top5 <- list("Adipocyte"  = c ("PLIN1", "LPL", "GPAM", "AQP7","ADIPOQ"),
                             "Venular VECs" = c("SELP", "NOSTRIN", "MYRIP", "GNA14", "TACR1"),
                             "Arteliolar VECs" = c("PODXL", "EFNB2", "NEBL", "SYT1", "THSD4"),
                             "vSMCs" = c("MYH11", "TAGLN", "SORBS2", "RCAN2", "LMOD1", "CARMN"),
                             "LEC" = c( "MMRN1", "PROX1", "PKHD1L1", "FLT4", "NRG3"),
                             "Pericytes" = c("PDGFRB","COL25A1", "POSTN", "STEAP4", "RGS5"),
                             "Nerve" = c("NRXN1", "XKR4","BMS1P14","CADM2","IL1RAPL2")
)

annotation_dotplot(so.stromal, stromal_markers_top5, "Stromal_markers", 15, 5)

muscle_markers_top5 <- list("Fast" = c("TNNT3", "MYH1", "TNNI2","MYBPC2", "ATP2A1"), 
                            "Trans" = c("COL22A1", "ADAMTSL1", "CPM", "SORBS2", "SPAG16"), 
                            "Slow"  = c("TNNC1", "TNNI1", "TNNT1", "ATP2A2", "LGR5"),
                            "Satellite" = c("PAX7", "CALCR", "CDH4", "CLCN5", "CTNND2"))
annotation_dotplot(so.muscle, muscle_markers_top5, "Muscle_markers", 12, 4)

#------------------------------------------------
# Figure S3D: FeaturePlots of key markers (broad)
#------------------------------------------------

broad_features <- c("COL1A2", "TTN", "PECAM1", "CD163", "NOTCH3", 
                    "ADIPOQ", "THEMIS", "MMRN1", "PAX7", "KIT", "FCRL1", 
                    "NRXN1", "JCHAIN")
do_FeaturePlot(so.achilles, reduction = "umap.scvi", features = broad_features,
               plot_cell_borders = TRUE, pt.size = 0.5, 
               individual.titles = broad_features,
               plot.title.face = "italic", ncol = 3,
               order = TRUE
               )

ggsave(paste0(directory, "/Supplementary/Figure_S2/Key_Featureplots.png"), width = 16, height = 24)




#------------------------------------------------
# Figure S4A: UMAP of Xenium broad annotation
# Figure S4B: Cell proportions by technology  
#------------------------------------------------

# Other script

#---------------------------------------------------
# Figure S4C: Bar chart of cell composition (broad)
#-----------------------------------------------------

# count the number of cells per celltype, microanatomy and patient
df <- so.achilles[[]] %>% select(broad_annotation, microanatomical_site, patient)
df <- data.table(df)
df <- df[, .(COUNT = .N), by = names(df)]

# calculate the proportion on a per-patient basis
list <- df %>% group_split(patient)

prop_fun <- function(df){
    df %>% mutate(proportion = COUNT/sum(COUNT)*100)
}

list <- lapply(list, prop_fun)

# find the mean proportion
df <- bind_rows(list)

df <- df %>% group_by(broad_annotation, microanatomical_site) %>% 
    dplyr::summarise(Mean = mean(proportion))

df$microanatomical_site <- str_replace(df$microanatomical_site, "Enth", "Enthesis")
df$microanatomical_site <- str_replace(df$microanatomical_site, "muscle", "Muscle")
df$microanatomical_site <- str_replace(df$microanatomical_site, "MB", "Midbody")
names(ma.cols) <- c("Enthesis", "Midbody", "MTJ", "Muscle")

celltypes <- c("Fibroblasts", 
               "Skeletal muscle cells", 
               "Vascular endothelial cells", 
               "Macrophages", 
               "Mural cells", 
               "Adipocytes", 
               "T cells", 
               "Lymphatic endothelial cells", 
               "Satellite cells", 
               "Granulocytes", 
               "B cells", 
               "Schwann cells", 
               "Plasma cells"
               )

df$broad_annotation <- factor(df$broad_annotation, levels = celltypes)
ggplot(df, aes(fill=microanatomical_site, y=Mean, x=broad_annotation)) + 
    geom_bar(position="fill", stat="identity")+
    theme_classic()+
    scale_fill_manual(values = ma.cols)+
    labs(x = "", y = "Proportion")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    labs(fill = "Microanatomical site")

ggsave(paste0(directory, "/Supplementary/Figure_S4/Cell_composition_barchart.png"), 
       width = 8, height = 5, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S4/Cell_composition_barchart.svg"), 
       width = 8, height = 5, bg = "white")
#---------------------------------------------------
# Figure S4D: Bar chart of cell composition (broad) in Xenium
#-----------------------------------------------------

# Other script

#---------------------------------------------------
# Figure S5A: Number of DE genes
#-----------------------------------------------------

# How many DE genes were in each cell type compared to Enthesis?

DE_Enth <- read.table("20250304_17-29_Pseudobulk.dir/DEseq2_results/Enth_results_summary.txt", 
                      header = TRUE, sep= "\t")
# pivot longer
DE_Enth.long <- DE_Enth %>% pivot_longer(3:5, values_to = "count", names_to = "type")
# merge column for cell annotation & comparison
DE_Enth.long <- DE_Enth.long %>% mutate (cell_annotation_comparison  = paste(cell_annotation_pseudobulk, comparison, sep = "_"))
DE_Enth.long

DE_Enth.long$broad_annotation_comparison <- factor(DE_Enth.long$cell_annotation_comparison, 
                                                   levels = c("Fibroblasts_MBvsEnth", "Fibroblasts_MTJvsEnth", "Fibroblasts_musclevsEnth", 
                                                              "Skeletalmuscle_MBvsEnth", "Skeletalmuscle_MTJvsEnth", "Skeletalmuscle_musclevsEnth",
                                                              "Adipocytes_MBvsEnth", "Adipocytes_MTJvsEnth", "Adipocytes_musclevsEnth",
                                                              "VEC_MBvsEnth", "VEC_MTJvsEnth", "VEC_musclevsEnth", 
                                                              "LEC_MBvsEnth", "LEC_MTJvsEnth", "LEC_musclevsEnth", 
                                                              "Mural_MBvsEnth", "Mural_MTJvsEnth", "Mural_musclevsEnth", 
                                                              "Granulocyte_MBvsEnth", "Granulocyte_MTJvsEnth", "Granulocyte_musclevsEnth", 
                                                              "Macrophages_MBvsEnth", "Macrophages_MTJvsEnth", "Macrophages_musclevsEnth", 
                                                              "Tcells_MBvsEnth", "Tcells_MTJvsEnth", "Tcells_musclevsEnth"
                                                   )
)

tol_muted <- c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", "#999933", "#AA4499", "#DDDDDD")
comparison.cols <-  c(tol_muted[2], tol_muted[4], tol_muted[6])
ggplot(DE_Enth.long, aes(x = broad_annotation_comparison, y = count))+
    geom_point(aes(colour = comparison, shape = type), size = 3)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    scale_colour_manual(values = c("MBvsEnth" = comparison.cols[1], "MTJvsEnth" = comparison.cols[2], "musclevsEnth"=comparison.cols[3]))+
    theme(axis.title.x = element_blank())+
    ylab("Number of differentially expressed genes")+
    theme(axis.text.y = element_text(size = 12, colour = "black"))+
    theme(axis.title.y = element_text(size = 12, colour = "black"))
ggsave(paste0(directory, "/Supplementary/Figure_S5/Number_DEgenes_vs_Enth.png"), width = 12, height = 6)
ggsave(paste0(directory, "/Supplementary/Figure_S5/Number_DEgenes_vs_Enth.svg"), width = 12, height = 6)

# How many DE genes were in each cell type compared to MTJ?

DE_MTJ <- read.table(paste0("20250304_17-29_Pseudobulk.dir/DEseq2_results/MTJ_results_summary.txt"), 
                     
                     header = TRUE, sep= "\t")

DE_MTJ.long <- DE_MTJ %>% pivot_longer(3:5, values_to = "count", names_to = "type")
# merge column for cell annotation & comparison
DE_MTJ.long <- DE_MTJ.long %>% mutate (cell_annotation_comparison = paste(cell_annotation_pseudobulk, comparison, sep = "_"))

DE_MTJ.long$cell_annotation_comparison <- factor(DE_MTJ.long$cell_annotation_comparison, 
                                                 levels = c("Fibroblasts_EnthvsMTJ", "Fibroblasts_MBvsMTJ", "Fibroblasts_musclevsMTJ",
                                                            "Skeletalmuscle_EnthvsMTJ", "Skeletalmuscle_MBvsMTJ", "Skeletalmuscle_musclevsMTJ", 
                                                            "Adipocytes_EnthvsMTJ", "Adipocytes_MBvsMTJ", "Adipocytes_musclevsMTJ", 
                                                            "VEC_EnthvsMTJ", "VEC_MBvsMTJ", "VEC_musclevsMTJ", 
                                                            "LEC_EnthvsMTJ", "LEC_MBvsMTJ", "LEC_musclevsMTJ", 
                                                            "Mural_EnthvsMTJ", "Mural_MBvsMTJ", "Mural_musclevsMTJ", 
                                                            "Granulocyte_EnthvsMTJ","Granulocyte_MBvsMTJ", "Granulocyte_musclevsMTJ", 
                                                            "Macrophages_EnthvsMTJ", "Macrophages_MBvsMTJ" , "Macrophages_musclevsMTJ", 
                                                            "Tcells_EnthvsMTJ", "Tcells_MBvsMTJ", "Tcells_musclevsMTJ"
                                                 )
)


ggplot(DE_MTJ.long, aes(x = cell_annotation_comparison, y = count))+
    geom_point(aes(colour = comparison, shape = type), size = 3)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    scale_colour_manual(values = c("EnthvsMTJ" = comparison.cols[1], "MBvsMTJ" = comparison.cols[2], "musclevsMTJ"=comparison.cols[3]))+
    theme(axis.title = element_blank())+
    ylab("Number of differentially expressed genes")+
    theme(axis.text.y = element_text(size = 12, colour = "black"))+
    theme(axis.title.y = element_text(size = 12, colour = "black"))
ggsave(paste0(directory, "/Supplementary/Figure_S5/Number_DEgenes_vs_MB.png"), width = 12, height = 6)
ggsave(paste0(directory, "/Supplementary/Figure_S5/Number_DEgenes_vs_MB.svg"), width = 12, height = 6)


df <- DE_Enth.long %>% group_by(cell_annotation_pseudobulk) %>% tally(count) 
df$cell_annotation_pseudobulk
df$cell_annotation <- c("Adipocytes", "Fibroblasts", "Granulocytes", "Lymphatic endothelial cells", 
                        "Macrophages", "Mural cells", "Skeletal muscle cells", 
                        "T cells", "Vascular endothelial cells")    
my_levels <- df %>% arrange(desc(n)) %>% pull(cell_annotation)
df$cell_annotation <- factor(df$cell_annotation, levels = my_levels)

ggplot(df, aes(x = cell_annotation, y = n, fill = cell_annotation))+
    geom_bar(stat = "identity")+
    theme_classic()+
    theme(axis.title = element_blank())+
    coord_flip()+
    scale_fill_manual(values = achilles.colours )+
    theme(axis.text.y = element_text(size = 20))+
    theme(axis.text.x = element_text(size = 14))+
    theme(legend.position="none")
ggsave(paste0(directory, "/Supplementary/Figure_S5/Number_DEgenes.png"), width = 10, height = 6)



#-----------------------------------------
# Figure S5D: Fibroblast annotation dotplot
#-----------------------------------------

annotation_dotplot <- function (so, genes){
    p <- DotPlot(so, features = genes, assay = "soupX") + 
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "#f7f3f2", mid = "#f2390a", high = "#6e1e0a", midpoint = 1)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16))+
        theme(axis.text.y = element_text(size = 16))+
        theme(axis.title.x = element_blank())+
        theme(axis.title.y = element_blank())+
        theme(strip.text = element_blank())
    
}

fibroblast_markers_top3 <- list("NEGR1" = c("NEGR1", "DCLK1", "BMP5", "FBN1"), 
                                "COMP" = c("COMP", "ITGA10", "KLHL29", "COL12A1"), 
                                "VCAN" = c("VCAN", "CACNB2", "ISM1"),
                                "MMP3" = c("MMP3", "FBLN1", "CILP"),
                                "COL15A1" = c("COL15A1", "GREB1L", "COL4A1", "ALDH1A2", "HSPG2", "APOD"),
                                "Chondrocyte" = c("ACAN", "COL2A1", "SDK2", "WWP2"),
                                "ITGA6" = c("TENM2", "ITGA6", "ABCA10", "CNTN4","SLC22A3"),
                                "THBS4" = c("THBS4", "NCAM1", "PIEZO2", "COL22A1"),
                                "PRG4" = c("PRG4", "CMKLR2", "ITGB8")) 
fibroblast_markers_top3 <- unlist(fibroblast_markers_top3, use.names = FALSE)

annotation_dotplot(so.fibroblast, fibroblast_markers_top3)
ggsave(paste0(directory, "/Supplementary/Figure_S5/Fibroblast_annotation_dotplot.png"), bg = "white", 
       width = 14, height = 5)

# change the order of the fibroblast names so they are grouped together

so.fibroblast$cell_annotation_pseudobulk <- factor(so.fibroblast$cell_annotation_pseudobulk, 
                                                   levels = c("NEGR1.VCAN.fb",
                                                              "NEGR1.COL15A1.fb", 
                                                              "NEGR1.ITGA6.fb", 
                                                              "COMP.MMP3.fb", 
                                                              "COMP.THBS4.fb", 
                                                              "Chondrocytes", 
                                                              "PRG4.fb"))
Idents(so.fibroblast) <- so.fibroblast$cell_annotation_pseudobulk




ggsave(paste0(directory, "/Supplementary/Figure_S5/Fibroblast_annotation_dotplot_2.svg"), bg = "white", 
       device = svg, width = 14, height = 5)


#----------------------------------------------------------------
# Figure S5E: UMAP of fibroblasts in Xenium
#----------------------------------------------------------------

# other scripts


#-----------------------------------------------
# Figure S5F: UMAP of two groups NEGR1+ and COMP+
#-----------------------------------------------


do_DimPlot(so.fibroblast, reduction = "umap", group.by = "fine_annotation_fibroblasts",
           colors.use = fibroblast.colours, 
           plot_cell_borders = TRUE, pt.size = 0.5, 
           legend.position = "right")
ggsave(paste0(directory, "/Supplementary/Figure_S4//UMAP_fibroblasts.png"), width = 10, height = 7, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S4//UMAP_fibroblasts.svg"), width = 10, height = 7, bg = "white")

#-----------------------------------------------
# Figure S5G: Broad fibroblast features
#-----------------------------------------------

broad_fibroblast_features <- c("NEGR1", "COMP")
do_FeaturePlot(so.fibroblast, reduction = "umap", features = broad_fibroblast_features,
               plot_cell_borders = TRUE, pt.size = 0.5, 
               individual.titles = broad_fibroblast_features,
               plot.title.face = "italic", ncol = 2,
               font.size = 20,
               order = TRUE
)

ggsave(paste0(directory, "/Supplementary/Figure_S5/NEGR1_COMP_Featureplots.png"), width = 10, height = 5)

#----------------------------------------------------------------
# Figure S5H: Feature plots of key genes (fibroblasts)
#----------------------------------------------------------------

fibroblast_features <- c("VCAN", "MMP3", "COL15A1", "ACAN", "ITGA6", "THBS4", "PRG4")
do_FeaturePlot(so.fibroblast, reduction = "umap", features = fibroblast_features,
               plot_cell_borders = TRUE, pt.size = 0.5, 
               individual.titles = fibroblast_features,
               plot.title.face = "italic", font.size = 20,
               ncol = 4,
               order = TRUE
)

ggsave(paste0(directory, "/Supplementary/Figure_S5/Fibroblast_Featureplots.png"), 
       width = 15, height = 8)

#---------------------------------------------------
# Figure S5I: Pathways enriched in fibroblast subsets
#---------------------------------------------------

# read in results
pathways.fb.subset <- read.table("20250305_13-28_Pathway_analysis_fibroblasts.dir/Pseudobulk/Combined_results.txt",
                                 sep = "\t", header = TRUE)
pathways.fb.subset$cluster <- factor(pathways.fb.subset$cluster, levels = unique(pathways.fb.subset$cluster))
# make IDs a factor
pathways <- as.data.frame(sort(table(pathways.fb.subset$ID)))
colnames(pathways) <- c("ID", "pathway_freq")
pathways.fb.subset <- pathways.fb.subset %>% 
    left_join(pathways) %>%
    group_by(pathway_freq) %>%
    arrange(cluster)

pathways.fb.subset$ID <- factor(pathways.fb.subset$ID, levels = rev(unique(pathways.fb.subset$ID)))    

# plot
ggplot(pathways.fb.subset, aes(x = cluster, y = ID, size = Count, colour = new_col)) +
    geom_point()+
    #scale_colour_viridis()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_blank())+
    labs(colour = "-log10(p.adj)")+
    scale_y_discrete(position = "right")
ggsave (paste0(directory, "/Supplementary/Figure_S5/Fibroblasts_subsets_pathways.png"), 
        width = 8, height = 10)
ggsave (paste0(directory, "/Supplementary/Figure_S5/Fibroblasts_subsets_pathways.svg"), 
        width = 8, height = 10)

#-----------------------------------------------
# Figure S5J: Bar chart of cell composition (fibroblasts) by technology
#-----------------------------------------------

# Xenium script 

#-----------------------------------------------
# Figure S5K: Bar chart of cell composition (fibroblasts)
#-----------------------------------------------

# count the number of cells per celltype, microanatomy and patient
df <- so.fibroblast[[]] %>% select(fine_annotation_fibroblasts, microanatomical_site, patient)
df <- data.table(df)
df <- df[, .(COUNT = .N), by = names(df)]

# calculate the proportion on a per-patient basis
list <- df %>% group_split(patient)

prop_fun <- function(df){
    df %>% mutate(proportion = COUNT/sum(COUNT)*100)
}

list <- lapply(list, prop_fun)

# find the mean proportion
df <- bind_rows(list)

df <- df %>% group_by(fine_annotation_fibroblasts, microanatomical_site) %>% 
    dplyr::summarise(Mean = mean(proportion))

df$microanatomical_site <- str_replace(df$microanatomical_site, "Enth", "Enthesis")
df$microanatomical_site <- str_replace(df$microanatomical_site, "muscle", "Muscle")
df$microanatomical_site <- str_replace(df$microanatomical_site, "MB", "Midbody")
names(ma.cols) <- c("Enthesis", "Midbody", "MTJ", "Muscle")

df$fine_annotation_fibroblasts <- factor(df$fine_annotation_fibroblasts, levels = c("COMPhi MMP3hi fibroblasts", 
                                                                                    "NEGR1hi VCANhi fibroblasts", 
                                                                                    "NEGR1hi COL15A1hi fibroblasts", 
                                                                                    "NEGR1hi ITGA6hi fibroblasts", 
                                                                                    "COMPhi THBS4hi fibroblasts", 
                                                                                    "Chondrocytes", 
                                                                                    "PRG4hi fibroblasts"))
ggplot(df, aes(fill=microanatomical_site, y=Mean, x=fine_annotation_fibroblasts)) + 
    geom_bar(position="fill", stat="identity")+
    theme_classic()+
    scale_fill_manual(values = ma.cols)+
    labs(x = "", y = "Proportion")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    labs(fill = "Microanatomical site")

ggsave(paste0(directory, "/Supplementary/Figure_S5/Cell_composition_barchart.png"), 
       width = 8, height = 5, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S5/Cell_composition_barchart.svg"), 
       width = 8, height = 5, bg = "white")

#-----------------------------------------------
# Figure S5L: Bar chart of cell composition (fibroblasts) in Xenium
#-----------------------------------------------

# Xenium script

#-----------------------------------------------
# Figure S5M: Examples of region-specific gene expression
#-----------------------------------------------


select_ma_genes <- c("GPC6", "MMP3", "PRG4", "TTN", "GREB1L")

ma.long <- data.frame(so.fibroblast$microanatomical_site)
ma.long$microanatomy <- str_replace(ma.long$so.fibroblast.microanatomical_site, "Enth", "Enthesis")
ma.long$microanatomy <- str_replace(ma.long$microanatomy, "muscle", "Muscle")
ma.long$microanatomy <- str_replace(ma.long$microanatomy, "MB", "Midbody")
names(ma.cols) <- c("Enthesis", "Midbody", "MTJ", "Muscle")
so.fibroblast <- AddMetaData(so.fibroblast, ma.long$microanatomy, col.name = "microanatomy")


# plot by microanatomy
VlnPlot(so.fibroblast, features = select_ma_genes, group.by = "microanatomy", 
        pt.size = 0, cols = ma.cols, log = FALSE, assay = "soupX", ncol = 5)+
    theme(strip.text.x = element_text(face = "italic"))

ggsave(paste0(directory, "/Supplementary/Figure_S5/ma_specific_genes/select_genes_by_ma.png"), width = 15, height =3)
ggsave(paste0(directory, "/Supplementary/Figure_S5/ma_specific_genes/select_genes_by_ma.svg"), width = 15, height =3, device = grDevices::svg)


# plot by celltype
VlnPlot(so.fibroblast, features = select_ma_genes, group.by = "cell_annotation_pseudobulk", 
        pt.size = 0, cols = fibroblast.colours.short, log = FALSE, assay = "soupX", ncol = 5)
ggsave(paste0(directory, "/Supplementary/Figure_S5/ma_specific_genes/select_genes_by_celltype.png"), width = 15, height =4)
ggsave(paste0(directory, "/Supplementary/Figure_S5/ma_specific_genes/select_genes_by_celltype.svg"), width = 15, height =4, device = grDevices::svg)


# -----------------------------------------------
# Figure S5N: Examples of region-specific gene expression
#-----------------------------------------------

# SCENIC script 

#-----------------------------------------------------
# Figure S8: Details of interactions with fibroblasts
#-----------------------------------------------------


liana_trunc <- liana_res %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

png(file=paste0(directory, "/Supplementary/Figure_S8/Liana_heatmap.png"), width=800, height=800)
heat_freq(liana_trunc)
dev.off()





##########################################
## Not currently included in manuscript ##
##########################################

# Violin plots for QC

SCpubr::do_ViolinPlot(so.achilles, 
                      features = c("decontX_contamination", "soupX_fraction"),
                      group.by = "orig.ident", 
                      ncol = 1, 
                      legend.position = "none")+
    theme(axis.title.x=element_blank())

ggsave(paste0(directory, "/Supplementary/Figure_S1/Vln_plots.png"), width = 12, height = 12)

do_FeaturePlot(so.achilles, reduction = "umap.scvi", features = c("decontX_contamination", "soupX_fraction"),
               plot_cell_borders = TRUE, pt.size = 0.5)
ggsave(paste0(directory, "/Supplementary/Figure_S1/UMAP_decontX_soupX.png"), width = 12, height = 7)

#---------------------
# MiloR beeswarm plots
#---------------------

# MiloR on fine annotation
da_results <- read.table("20250226_10-11_MiloR_Cell_Proportions.dir/Results.dir/Achilles-milo-DA-results.txt", 
                         sep = "\t", header = TRUE)
da_results <- da_results %>% filter(fine_annotation != "Mixed")
da_results$fine_annotation <- factor(da_results$fine_annotation, 
                                      levels = c("Satellite cells",
                                                 "Slow-twitch skeletal muscle cells", 
                                                 "Fast-twitch skeletal muscle cells", 
                                                 "Transitional skeletal muscle cells",
                                                 "Pericytes", 
                                                 "NEGR1hi COL15A1hi fibroblasts", 
                                                 "COMPhi THBS4hi fibroblasts", 
                                                 "MERTKhi LYVE1hi macrophages", 
                                                 "Adipocytes", 
                                                 "Lymphatic endothelial cells", 
                                                 "Plasma cells", 
                                                 "Granulocytes", 
                                                 "CLEC10A DCs",                        
                                                 "Monocytes",
                                                 "B cells", 
                                                 "NEGR1hi ITGA6hi fibroblasts", 
                                                 "PRG4hi fibroblasts", 
                                                 "Schwann cells",               
                                                 "T cells",
                                                 "vSMC", 
                                                 "Venular VEC", 
                                                 "Arteriolar VEC", 
                                                 "COMPhi MMP3hi fibroblasts",
                                                 "NEGR1hi VCANhi fibroblasts",
                                                 "Chondrocytes"  ))
plotDAbeeswarm(da_results, group.by = "fine_annotation")+
    theme(axis.title.y=element_blank())+
    scale_fill_gradient2(low = "#00466B", mid = "#FFFFFF", high = "#A91400", midpoint = 0, name = "logFC")
#ggsave(paste0(directory, "/Supplementary/Beeswarm_plots/Beeswarm_plot_fine.png"), width = 10, height = 12, bg = "white")
#ggsave(paste0(directory, "/Supplementary/Beeswarm_plots/Beeswarm_plot_fine.svg"), width = 10, height = 10, bg = "white")



#---------------------------------------------------------
# Expression of genes related to Achilles tendinopathy GWAS
#---------------------------------------------------------

annotation_dotplot <- function (so, genes){
    p <- DotPlot(so, features = genes, assay = "soupX") + 
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "#f7f3f2", mid = "#f2390a", high = "#6e1e0a", midpoint = 1)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16))+
        theme(axis.text.y = element_text(size = 16))+
        theme(axis.title.x = element_blank())+
        theme(axis.title.y = element_blank())+
        theme(strip.text = element_blank())
    
}

AT_GWAS_loci <- c("COL15A1", "CASP8", "CDCP1", "MPP7", "SOX21", 
                  "ADAMTS14", "BMP4", "FBN2", "MIR608", "MMP3", "TIMP2", "TNC")
p <- annotation_dotplot(so.fibroblast, AT_GWAS_loci)
p

T2D_GWAS_loci <- c("MTOR","TTN","VEGFA","PDGFC","ADAMTS9","PPARG",
                   "FGFR2","FGFR4","ACVR1C","ITGA1","LEPR","IRS1","IGF2","IGF2BP2")

p <- annotation_dotplot(so.fibroblast, T2D_GWAS_loci)
p

T2D_loci <- read.table("Files.dir/T2D_loci")

# top100 markers per cell type
broad_markers_top100 <- broad_markers %>%
    group_by(cluster) %>% 
    slice_max(n=100, order_by = abs(avg_log2FC)) %>% 
    pull(gene)

my_genes <- broad_markers_top100[which(broad_markers_top100 %in% T2D_loci$V1)]
head(broad_markers_top100)
head(T2D_loci$V1)
"NEGR1" %in% T2D_loci$V1
"NEGR1" %in% broad_markers_top100
"NEGR1" %in% my_genes
sort(my_genes)

# which cell types are they from?
broad_markers[which(broad_markers$gene %in% my_genes),] 

# Plot genes from fibroblasts
T2D_loci_fibroblasts <- broad_markers[which(broad_markers$gene %in% my_genes),] %>% 
    filter(cluster == "Fibroblasts") %>% pull(gene)

p <- annotation_dotplot(so.fibroblast, T2D_loci_fibroblasts)
p

p <- annotation_dotplot(so.achilles, T2D_loci_fibroblasts)
p


#------------------------------------------------
# Not used: Cluster plots of DE gene expression
#------------------------------------------------

celltypes <- c("Adipocytes", 
               "Fibroblasts",
               "Macrophages",
               "Mural",
               "Skeletalmuscle",
               "VEC")

clusterplot_list <- list()
for (cell in celltypes){
    clusterplot_list[[cell]] <- readRDS(paste0("20250304_17-29_Pseudobulk.dir/Cluster_results/", cell, "_plot.rds"))
}

clusterplot_list[["VEC"]] +
    theme_classic()+
    theme(legend.position = "none")+
    scale_colour_manual(values = "#00466B")+
    ggtitle("VEC")
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/VEC_clusters.png"), width = 8, height = 4)
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/VEC_clusters.svg"), width = 8, height = 4)

clusterplot_list[["Adipocytes"]] +
    theme_classic()+
    theme(legend.position = "none")+
    scale_colour_manual(values = "#00466B")+
    ggtitle("Adipocytes")
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/Adipocytes_clusters.png"), width = 10, height = 4)
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/Adipocytes_clusters.svg"), width = 10, height = 4)

clusterplot_list[["Macrophages"]] +
    theme_classic()+
    theme(legend.position = "none")+
    scale_colour_manual(values = "#00466B")+
    ggtitle("Macrophages")
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/Macrophages_clusters.png"), width = 10, height = 4)
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/Macrophages_clusters.svg"), width = 10, height = 4)

clusterplot_list[["Mural"]] +
    theme_classic()+
    theme(legend.position = "none")+
    scale_colour_manual(values = "#00466B")+
    ggtitle("Mural")
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/Mural_clusters.png"), width = 10, height = 4)
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/Mural_clusters.svg"), width = 10, height = 4)

clusterplot_list[["Skeletalmuscle"]] +
    theme_classic()+
    theme(legend.position = "none")+
    scale_colour_manual(values = "#00466B")+
    ggtitle("Skeletal muscle")
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/Skeletalmuscle_clusters.png"), width = 5, height = 4)
#ggsave(paste0(directory, "/Supplementary/Figure_S9/DEseq2/Skeletalmuscle_clusters.svg"), width = 5, height = 4)


#-----------------------------------------
# Expression of some interaction genes
#------------------------------------------

# Collagen-Integrin
VlnPlot(so.fibroblast, features = c("COL1A2", "FN1", "ITGA1"), group.by = "microanatomical_site", 
        pt.size = 0, cols = ma.cols, log = FALSE, assay = "soupX")
VlnPlot(so.fibroblast, features = c("COL1A2","FN1", "ITGA1", "ITGA5", "ITGAV", "ITGB1"), group.by = "fine_annotation_fibroblasts", 
        pt.size = 0, log = FALSE, assay = "soupX")
#ggsave(paste0(directory, "/Supplementary/Figure_S9/Collagens-integrins.png"), width = 15, height =10)

# TGFB
VlnPlot(so.fibroblast, features = c("TGFB2", "TGFBR3"), group.by = "microanatomical_site", 
        pt.size = 0, cols = ma.cols, log = FALSE, assay = "soupX")
VlnPlot(so.fibroblast, features = c("TGFB2", "TGFBR3"), group.by = "fine_annotation_fibroblasts", 
        pt.size = 0, log = FALSE, assay = "soupX")



# plot by celltype
VlnPlot(so.fibroblast, features = select_ma_genes, group.by = "cell_annotation_pseudobulk", 
        pt.size = 0, cols = fibroblast.colours.short, log = FALSE, assay = "soupX", ncol = 5)
#ggsave(paste0(directory, "/Figure_4/ma_specific_genes/select_genes_by_celltype.png"), width = 15, height =4)
#ggsave(paste0(directory, "/Figure_4/ma_specific_genes/select_genes_by_celltype.svg"), width = 15, height =4, device = grDevices::svg)


#---------------------------------------------------
# Figure: Bar chart of cell composition (muscle)
#-----------------------------------------------------

df <- so.muscle[[]] %>% select(fine_annotation_muscle, microanatomical_site)
df <- data.table(df)
df <- df[, .(COUNT = .N), by = names(df)]
df$microanatomical_site <- str_replace(df$microanatomical_site, "Enth", "Enthesis")
df$microanatomical_site <- str_replace(df$microanatomical_site, "muscle", "Muscle")
df$microanatomical_site <- str_replace(df$microanatomical_site, "MB", "Midbody")
names(ma.cols) <- c("Enthesis", "Midbody", "MTJ", "Muscle")


ggplot(df, aes(fill=microanatomical_site, y=COUNT, x=fine_annotation_muscle)) + 
    geom_bar(position="fill", stat="identity")+
    theme_classic()+
    scale_fill_manual(values = ma.cols)+
    labs(x = "", y = "Proportion")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    labs(fill = "Microanatomical site")

ggsave(paste0(directory, "/Supplementary/Figure_S7/Cell_composition_barchart_muscle.png"), 
       width = 8, height = 5, bg = "white")
ggsave(paste0(directory, "/Supplementary/Figure_S7/Cell_composition_barchart_muscle.svg"), 
       width = 8, height = 5, bg = "white")

#-------------------------------------------------------
# Record session Info
#-------------------------------------------------------


sink(paste0(directory, "/sessionInfo.txt"))
sessionInfo()
sink()


