# Hierarchical-Maker-Genes-Selection-for-scRNA-seq-Data
This is an executable R Code for paper Hierarchical Marker Genes Selection for scRNA-seq Data with the data used included.

The data used in this example is the raw PBMC3k dataset, located in filtered_gene_bc_matrix folder. We apply a standard preprocessing process for PBMC3k dataset in our implements

To construct the hierarchy in our proposed method, run pbmc_3k_hierarchy.R

To generate the compact heatmap based on the hierarchy, run Generate_compact_heatmap.R(This file should be run after constructing the hierarchy, which means it should be run after running pbmc_3k_hierarchy.R or with an predefined hierarchy organized in a proper data structure)

All the novel functions applied is included in Hierachy_utilities.R 

## SessionInfo():

R version 4.2.1 (2022-06-23 ucrt)\
Running under: Windows 10 x64 (build 19044)

### Attached base packages:

grid\
stats\
graphics\
grDevices\
utils\
datasets\
methods\
base

### other attached packages:

ggplot2_3.4.0\
ComplexHeatmap_2.12.1\
stringr_1.5.0\
reshape_0.8.9\
patchwork_1.1.2\
magick_2.7.5\
tidyseurat_0.7.2\
ttservice_0.3.8\
SeuratObject_4.1.3\
Seurat_4.3.0\
dplyr_1.1.0
