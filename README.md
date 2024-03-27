# Hierarchical-Maker-Genes-Selection-for-scRNA-seq-Data
This is an executable R Code for paper Hierarchical Marker Genes Selection for scRNA-seq Data with the data used included.

The data used in this example is the raw PBMC3k dataset, located in filtered_gene_bc_matrix folder. We apply a standard preprocessing process for PBMC3k dataset in our implements

To construct the hierarchy in our proposed method, run pbmc_3k_hierarchy.R

To generate the compact heatmap based on the hierarchy, run Generate_compact_heatmap.R(This file should be run after constructing the hierarchy, which means it should be run after running pbmc_3k_hierarchy.R or with an predefined hierarchy organized in a proper data structure)

All the novel functions applied is included in Hierachy_utilities.R 

Users can modify the number of marker genes selected to calculate the score(there are three places that need this parameter:function MergeClusters and BuildHierarchicalMap in Hierachy_utilities.R and after the FindAllMarkers function in pbmc_3k_hierarchy.R)

For PBMC3k dataset, we use the top10 marker genes for each lineage to generate the hierarchy and for PBMC control/stimulate we use top50 in terms of the larger number of genes.

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

Our algorithm are also compatible with the latest version Seurat v5. Users can choose to load Seurat v5 while taking care of the version of other dependencies for the following visualization and analysis. 

===================================================================================================================================================
# Evaluation
The evaluation is implemented on Python Jupyter notebook(Evaluation.ipynb), including data loading/preprocessing, classification and UMAP visualization, with detailed instructions included.

For the evaluation, we start from loading the raw data, which means the whole process is fully independent from the algorithm we construct the hierarchy implemented by R. This allows the users to customize their own way to process the data for evaluation in terms of different goals they want to achieve. 

To do the evaluation:

1. run pbmc_3k_hierarchy.R and get the hierarchy(the hierarchy structure and corresponding marker genes are saved as .csv files)
2. Open the Evaluation.ipynb and load the gene expression data, hierarchy structure and corresponding marker genes files we get from step 1(One can also manually input the hierarchy and the corresponded marker genes for each lineage into the python notebook, example codes are included)
4. follow the instructions in the notebook and change the parameters to do the evaluation.

### Python packages needed:

scanpy~=1.7.2\
scipy~=1.5.4\
pandas~=1.1.5\
numpy~=1.19.5\
sklearn~=0.0.post1\
scikit-learn~=0.24.1\
ismember~=1.0.1\
umap-learn~=0.5.1
