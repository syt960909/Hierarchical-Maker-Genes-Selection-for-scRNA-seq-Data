library(dplyr)
library(Seurat)
library(tidyseurat)

pbmc.data <- Read10X(data.dir = "Y:/qiu-lab/mapping_to_atlas/datasets/PBMC_10x_3K/filtered_gene_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

`%>%` = magrittr::`%>%`
pbmc_<-pbmc %>%  arrange(pbmc$seurat_clusters)
image(pbmc_[top10$gene])


pbmc_<-pbmc
#new.cluster.ids <- c("Naive CD4 T","CD14+ Mono","Memory CD4 T","B","CD8 T","FCGR3A+ Mono","NK","DC","Platelet") 
new.cluster.ids <- c("Naive CD4 T","CD14+ Mono","Memory CD4 T","B","CD8 T","FCGR3A+ Mono","NK","DC","Platelet") 
names(new.cluster.ids) <- levels(pbmc_)
pbmc_ <- RenameIdents(pbmc_, new.cluster.ids)
pbmc_ <- NormalizeData(pbmc_, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_ <- FindVariableFeatures(pbmc_, selection.method = "vst", nfeatures = 2000)
pbmc_ <- RunPCA(pbmc_, features = VariableFeatures(object = pbmc_))

# cluster re-assignment occurs, which re-assigns clustering in my_levels (assuming you have 12 clusters in total)
my_levels <- c('Naive CD4 T', 'Memory CD4 T', 'CD8 T', 'B', 'NK', 'CD14+ Mono', 'FCGR3A+ Mono', 'DC', 'Platelet')

# Re-level object@ident
pbmc_@active.ident <- factor(x = pbmc_@active.ident, levels = my_levels)
pbmc_ <- ScaleData(pbmc_, features = all.genes)###??????
pbmc.markers <- FindAllMarkers(pbmc_, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc_, features = top10$gene)


new.cluster.ids <- c("CD14 Mono","pDC","CD4 Memory T;CD4 Naive T","T activated","CD4 Memory T;CD4 Naive T","CD8 T;NK ","Mk","B activated;B","B activated;B","DC", "CD16 Mono","CD8 T;NK ", "Eryth")
names(new.cluster.ids) <- levels(PBMC_control)
PBMC_control <- RenameIdents(PBMC_control, new.cluster.ids)

names(new.cluster.ids) <- levels(PBMC_stim)
PBMC_stim <- RenameIdents(PBMC_stim, new.cluster.ids)
my_levels <- c("CD14 Mono;CD16 Mono;pDC;DC", "CD4 Memory T;T activated;CD4 Naive T;CD8 T;B activated;B;NK ","Mk","Eryth")
PBMC_control@active.ident <- factor(x = PBMC_control@active.ident, levels = my_levels)


pbmc_sub1 <- subset(x = pbmc, idents = c("Naive CD4 T", "Memory CD4 T"))

#new.cluster.ids<-c("Memory CD4 T;Naive CD4 T","CD14+ Mono;FCGR3A+ Mono;DC","Memory CD4 T;Naive CD4 T","B","NK;CD8 T","CD14+ Mono;FCGR3A+ Mono;DC","NK;CD8 T","CD14+ Mono;FCGR3A+ Mono;DC","Mk")
new.cluster.ids<-c("Naive CD4 T;Memory CD4 T", "Naive CD4 T;Memory CD4 T", "CD8 T")

names(new.cluster.ids) <- levels(pbmc_sub1)

pbmc_sub1 <- RenameIdents(pbmc_sub1, new.cluster.ids) 

pbmc_sub1 <- ScaleData(pbmc_sub1, features = all.genes)###??????
pbmc.markers <- FindAllMarkers(pbmc_sub1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc_sub1, features = top10$gene)
#GetHeatmapMatrix: to get the matrix used for plotting marker gene heatmap
data_<-GetHeatmapMatrix(pbmc,features = top10$gene)

##TODO:maybe start from average of each cell type and a list of genes, to calculate whether a heatmap looks good
score<-HeatmapScore(pbmc,data_,top10,levels(pbmc))
pbmc$best_score<-score
pbmc_3k_hie_list<-BuildHierarchicalMap(pbmc)


pbmc.markers <- FindAllMarkers(pbmc_3k_hie_list[[1]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc_3k_hie_list[[1]], features = top10$gene)