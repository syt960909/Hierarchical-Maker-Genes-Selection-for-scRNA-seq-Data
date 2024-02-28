library(dplyr)
library(Seurat)
library(tidyseurat)
source('Hierachy_utilities.R')  ###Contains the functions for building hierarchy
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
###Create heatmap for the current lineage
DoHeatmap(pbmc, features = top10$gene)
#GetHeatmapMatrix: to get the matrix used for plotting marker gene heatmap
data_<-GetHeatmapMatrix(pbmc,features = top10$gene)

###Get a initial heatmap score corresponding to a flat data structure
score<-HeatmapScore(pbmc,data_,top10,levels(pbmc))
pbmc$best_score<-score
###Build the hierarchy
pbmc_3k_hie_list<-BuildHierarchicalMap(pbmc)



#####################################################################################################
###Save results
###save hierarchical structure and marker genes corresponding to each hierarchical level as output files

###Hierarchical structure
hierarchy_list = list()
i = 2
while(i <= length(pbmc_3k_hie_list)){
  tmp_list = matrix(nrow = 0,ncol = length(levels(pbmc_3k_hie_list[[i]])))
  for(j in seq_along(levels(pbmc_3k_hie_list[[i]]))){
    print(j)
    tmp_list[j] = levels(pbmc_3k_hie_list[[i]])[j]
  }
  print(tmp_list)
  hierarchy_list = c(hierarchy_list,list(tmp_list))
  i = i + 1
}

hierarchy_gene_list = hierarchy_list
i = 2
while(i <= length(pbmc_3k_hie_list)){
  PBMC.markers <- FindAllMarkers(pbmc_3k_hie_list[[i]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top10 <- PBMC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  for(j in seq_along(hierarchy_list[[i-1]])){
    hierarchy_gene_list[[i-1]][[j]] = paste(top10$gene[top10$cluster==levels(pbmc_3k_hie_list[[i]])[[j]]],collapse = ";")
  }
  i = i+1
}

for(i in seq_along(hierarchy_list)){
  for(j in seq_along(hierarchy_list[[i]])){
    if(length(str_split_1(hierarchy_list[[i]][[j]],";"))==2){
      hierarchy_list = c(hierarchy_list,str_split(hierarchy_list[[i]][[j]],";"))
      hierarchy_gene_list = c(hierarchy_gene_list,str_split(hierarchy_list[[i]][[j]],";"))
      PBMC_sub <- subset(x = pbmc, idents = str_split_1(hierarchy_list[[i]][[j]],";"))
      all.genes <- rownames(PBMC_sub)
      PBMC_sub <- ScaleData(PBMC_sub, features = all.genes)
      PBMC.markers <- FindAllMarkers(PBMC_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      top10 <- PBMC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
      
      hierarchy_gene_list[[length(hierarchy_gene_list)]][[1]]=paste(top10$gene[top10$cluster==levels(PBMC_sub)[[1]]],collapse = ";")
      hierarchy_gene_list[[length(hierarchy_gene_list)]][[2]]=paste(top10$gene[top10$cluster==levels(PBMC_sub)[[2]]],collapse = ";")
      
      
    }
  }
}


###Please Modify the path to save the results
write.table(as.data.frame(hierarchy_list),file="your_own_path\\hierarchy_list.csv", quote=F,sep=",",row.names=F,col.names = F)
write.table(as.data.frame(hierarchy_gene_list),file="your_own_path\\hierarchy_gene_list.csv", quote=F,sep=",",row.names=F,col.names = F)
