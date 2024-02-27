library(dplyr)
library(Seurat)
library(tidyseurat)
source('Hierachy_utilities.R') ###Contains the functions for building hierarchy

#####################################################################################################
###This demo script provides an alternative faster implementation
###For the original implementation, we do wilcoxon test for every possible combination, which is time consuming
###Here we provide a faster implementation based on Seurat BuildClusterTree function, which constructs a phylogenetic tree relating the 'average' cell from each identity class.
###Tree is estimated based on a distance matrix constructed in either gene expression space or PCA space(default).
###Through the BuildClusterTree, we can get a vector called "height", which is a set of n-1 real values (non-decreasing for ultrametric trees). 
###The clustering height: that is, the value of the criterion associated with the clustering method for the particular agglomeration.
###From bottom to top, we set each height value as threshold to cut the tree into groups so that we get the hierarchical lineages
###We calculate the heatmap score corresponding to lineages for each level, and stop the process when the score is not getting higher
###We keep the lineages under the stop height as the final hierarchical tree. 



###Load data
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrix")


#####################################################################################################
###Data preprocessing with Seurat
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#####################################################################################################
###Data Clustering&annotation(if the cell type annotation is provided, we can skip this step)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

#####################################################################################################
###Using Seurat FindAllMarkers function to find marker genes in a one-vs-all fashion 
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
###Select top10 genes for each cluster
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(pbmc, features = top10$gene)
#GetHeatmapMatrix: to get the matrix used for plotting marker gene heatmap
data_<-GetHeatmapMatrix(pbmc,features = top10$gene)

###Get a initial heatmap score corresponding to a flat data structure
score<-HeatmapScore(pbmc,data_,top10,levels(pbmc))
pbmc$best_score<-score
pbmc_3k_hie_list_fast<-BuildHierarchicalMap_predefine(pbmc)
