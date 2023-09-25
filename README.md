# Hierarchical-Maker-Genes-Selection-for-scRNA-seq-Data
## This is an executable R Code for paper Hierarchical Marker Genes Selection for scRNA-seq Data with the data used included.
## The data used in this example is the raw PBMC3k dataset, located in filtered_gene_bc_matrix folder. We apply a standard preprocessing process for PBMC3k dataset in our implements/
## To construct the hierarchy in our proposed method, run pbmc_3k_hierarchy.R
## To generate the compact heatmap based on the hierarchy, run Generate_compact_heatmap.R(This file should be run after constructing the hierarchy, which means it should be run after running pbmc_3k_hierarchy.R or with an predefined hierarchy organized in a proper data structure)
## All the novel functions applied is included in Hierachy_utilities.R 
