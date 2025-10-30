# Fibroblast specialisation across microanatomy in a single-cell atlas of healthy human Achilles tendon
Carla J. Cohen,  Jolet Y. Mimpen,  Alina Kurjan, Claudia Paul, Shreeya Sharma, Lorenzo Ramos-Mucci, Chinemerem T. Ikwanusi, Ali Cenk Aksu, Tracy Boakye Serebour, Marina Nikolic, Kevin Rue-Albrecht, Christopher Gibbons, Duncan Whitwell, Tom Cosker, Steven Gwilym, Ather Siddiqi, Raja Bhaskara Rajasekaran, Harriet Branford-White, Adam P. Cribbs, Philippa A. Hulley, David Sims, Mathew J. Baldwin, Sarah J. B. Snelling 
  
[BioRxiv preprint](https://doi.org/10.1101/2025.10.17.683019)

Code by [Carla Cohen](https://github.com/carlacohen), [Jolet Mimpen](https://orcid.org/0000-0003-4464-242X) and [Alina Kurjan](https://orcid.org/0000-0002-4503-7865)

The following scripts were used to analyse single nucleus RNA-seq data and spatial transcriptomics data from human healthy Achilles tendon at three microanatomical sites, and adjoining muscle.  Scripts can be run using conda environments; the yaml files to create these are in the conda directory. Parameters for each script are set in corresponding yaml files in the yml directory.  


### QC, filtering and integration of single nuc RNA-seq data  
|  | Script | Purpose |
| ----- | ------ | ------- |
| 1 | QC-filter-achilles-combined.Rmd | Initial QC and filtering on a per sample basis |
| 2 |  Ambient-doublet-achilles-combined.Rmd | Remove ambient RNA with decontX and calculate doublets | 
| 3 |  SoupX-achilles-combined.Rmd| SoupX on individual samples using automated method | 
| 4 |  SoupX-manual-achilles-combined.Rmd | SoupX on individual samples using manual method | 
| 5 |  Doublet-achilles-combined.Rmd| Identify and remove doublets | 
| 6 | Integration-achilles-combined.Rmd | Harmony integration |
| 7 | Achilles_scVI_autotune.ipynb | Set hyperparameters for scVI |
| 8 | Achilles_scVI.ipynb  | Integration with scVI |
| 9 |  Integration-achilles-scVI-h5ad-visualisation.Rmd | Visualise outcome of scVI integration |


### Annotation of single nuc RNA-seq data

|  | Script | Purpose |
| ----- | ------ | ------- |
| 1 | Annotation-achilles-combined.Rmd | Broad annotation | 
| 2 | Annotation-update-achilles-combined-soupX_scVI_louvain.Rmd  | Update broad annotation |
| 3 | scVI_on_fibroblast_subset.ipynb | Subset and recluster fibroblasts with scVI | 
| 4 | scVI_on_immune_subset.ipynb | Subset and recluster immune with scVI | 
| 5 | scVI_on_stromal_subset.ipynb | Subset and recluster stromal with scVI | 
| 6 | Subset_recluster_Muscle.Rmd | Subset and recluster muscle subset with Harmony | 
| 7 | Achilles_Adipocyte.VEC.LEC.Mural_Fine-annotation_scVI_louvain_0.5.ambient.Rmd | Annotation of stromal subset |
| 8 | Achilles_Fibroblast_Fine-annotation.Rmd | Annotation of fibroblast subset | 
| 9 | Achilles_Immune_Fine-annotation_scVI_louvain_1.6_soupX.Rmd | Fine annotation of immune subset| 
| 10 | Achilles_Muscle_Fine-annotation_harmony_louvain_0.1_soupX.Rmd | Fine annotation of muscle subset | 
| 11 | Achilles_Stromal_Fine-annotation_scVI_louvain_0.9.ambient.Rmd | Fine annotation of stromal subset | 
| 12 | Achilles_final_annotation.Rmd| Bring together all fine annotations| 


### Downstream analysis of single nuc RNA-seq data
|  | Script | Purpose |
| ----- | ------ | ------- |
| 1 | MiloR_Cell_Proportions-Achilles.Rmd | MiloR analysis of cell proportions | 
| 2 | Pseudobulk_DEseq2_LRT_Achilles.Rmd | Pseudobulk and DEseq2 of broad annotation | 
| 3 | Pseudobulk_DEseq2_LRT_fibroblasts.Rmd | Pseudobulk and DEseq2 of fibroblast subset | 
| 4 | Pathway_analysis_Achilles.Rmd| Pathway analysis on broad annotation pseudobulk genes | 
| 5 | Fibroblasts_pathway_analysis.Rmd | Pathway analysis on fibroblast subset | 
| 6 | Liana-achilles-fine_annotation.Rmd | Calculate cell-cell interactions on fine annotation and circle plots | 
| 7 | Liana.R |Calculate cell-cell interaction scores on more databases | 
| 8 | Adule Ach SCENIC.ipynb | SCENIC analysis of fibroblast subset |

### Xenium spatial transcriptomics analysis  
The following scripts were used to analyse Xenium spatial transcriptomics data (in the Xenium directory)  
  
|  | Script | Purpose |
| ----- | ------ | ------- |
| 1 | Xenium_panel_designer.R | Prepare objects for panel design 
| 2 | Xenium_custom_panel.Rmd | Exploratory script to design Xenium custom panel | 
| 3 | Xenium_markers.Rmd | Expression of proposed custom panel on single cell datasets |
| 4 | Xenium_QC_metrics.Rmd | Plot QC metrics from the Xenium analyser | 
| 5 | Xenium_Achilles.Rmd | QC and filtering for each sample individually |
| 6 | Xenium_integration.Rmd | Integration of four samples |
| 7 | Cluster_annotation.Rmd| Annotation by clustering rather than label transfer (not used in manuscript) |
| 8 | Fibroblast_annotation.Rmd | Annotation of fibroblasts by label transfer |
| 9 | Immune_annotation.Rmd | Annotation of immune cells by label transfer |
| 10 | Muscle_annotation.Rmd | Annotation of muscle cells by label transfer |
| 11 | Stromal_annotation.Rmd | Annotation of stromal cells by label transfer | 
| 12 | Xenium_visualisation.Rmd | Harmonisation of annotations and visualisation | 
| 13 | Achilles_paper_figures_Xenium.R | Create figures for the manuscript |


### Convert objects between formats 
|  | Script | Purpose |
| ----- | ------ | ------- |
| 1 | Convert_h5ad_seurat.Rmd| Convert h5ad objects to seurat | 
| 2 | Convert_sce_h5ad-Achilles.R | Convert a list of sce objects to h5ad | 
| 3 | Convert_sce_h5ad.R | Convert a single sce object to h5ad | 
| 4 | Convert_seurat_h5ad.R | Convert seurat to h5ad format | 

### Other useful scripts  

General scripts to help with analysis but not essential for final outcome   

|  | Script | Purpose |
| ----- | ------ | ------- |
| 1 | Achilles_paper_figures.R | Generate paper figures | 
| 2 | Achilles_annotation_template.Rmd | Lists of annotation markers |
| 3 | Achilles_subset_annotations.Rmd | Plot all annotation markers on any subset | 
| 4 | Cluster-achilles-combined.Rmd| Clustering of individual samples for QC | 
| 5 |  Cluster-leiden-Achilles-combined.Rmd | Cluster with leiden after integration to compare with louvain (not used in manuscript) | 
| 6 | Compare_Harmony_scVI.Rmd | Additional QC for choosing best integration method | 
| 7 | Integration-QC-achilles.Rmd | Additional integration QC | 
| 8 | Liana-plotting.Rmd | Plots of cell-cell interactions | 
| 9 | Matrisome_Achilles_fibroblasts.Rmd | Analysis of matrisome gene expression (not used in manuscript) | 
| 10 | Subset_recluster.Rmd | Subset and recluster any group of cells | 
| 11 | Cluster_heterogeneity.R| Assess cluster heterogeneity using ROGUE | 



