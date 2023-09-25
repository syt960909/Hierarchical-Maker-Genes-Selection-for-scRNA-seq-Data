library(ComplexHeatmap)
if(!any(grepl("package:ggplot2", search()))) {
  message("loading ggplot2 for shinystan visualizations")
  library(ggplot2, quiet =TRUE)
}
library(ggplot2)
library(patchwork)
library(reshape)
library(magick)
library(stringr)
source('Hierachy_utilities.R')
my_levels<-character(0)
i = length(pbmc_3k_hie_list)
while (i>1){
  for (ele in levels(pbmc_3k_hie_list[[i]])){
    for (ele_sub in str_split(ele,';')[[1]]){
      
      if (!(ele_sub %in% my_levels)){
        my_levels<-c(my_levels,ele_sub)
      }
      
    }
    
  }
  i=i-1
}
###########re-order my_level##########
#check the hierarchy of first celltype
branch_list<-list()

for(my_level_ele in my_levels){
  i = length(pbmc_3k_hie_list)
  route<-character(0)
  while (i>1){
    
    for (ele in levels(pbmc_3k_hie_list[[i]])){
      
      if (my_level_ele %in% str_split(ele,';')[[1]]){
        print(str_split(ele,';')[[1]])
        route<-c(route,i)
        }
    }
   i=i-1
  }
  branch_list<-append(branch_list,list(route))
}

over_lap<-list()
for (i in seq_along(1:length(branch_list)-1)){
  print(i)
  over_lap<-append(over_lap,as.numeric(length(intersect(branch_list[[1]],branch_list[[i+1]]))))
  print(intersect(branch_list[[1]],branch_list[[i+1]]))
}
names(over_lap)<-my_levels[2:length(my_levels)]
over_lap_s<-over_lap[order(unlist(over_lap),decreasing=TRUE)]
my_levels[2:length(my_levels)]<-names(over_lap_s)
for(my_level_ele in my_levels){
  print(my_level_ele)
}


# Re-level object@ident
pbmc_3k_hie_list[[1]]@active.ident <- factor(x = pbmc_3k_hie_list[[1]]@active.ident, levels = my_levels)
pbmc.markers <- FindAllMarkers(pbmc_3k_hie_list[[1]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
a<-DoHeatmap_matrix(pbmc_3k_hie_list[[1]],features = top10$gene)
cell_order<-a$cell.order
col_anno <- HeatmapAnnotation(group = a$group.by,name = "ann") 
###########generate all sub-heatmaps#############################
ht_list = col_anno
i = 1
while (i<length(pbmc_3k_hie_list)){
  pbmc.markers <- FindAllMarkers(pbmc_3k_hie_list[[i+1]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  a<-DoHeatmap_matrix(pbmc_3k_hie_list[[i+1]],features = top10$gene)
  data_tmp <- data.frame(matrix(0,nrow = 2638, ncol = ncol(a$data)))
  colnames(data_tmp)<-colnames(a$data)
  rownames(data_tmp)<-colnames(pbmc@assays$RNA@counts)

  for (bar_code in colnames(pbmc@assays$RNA@counts)) {
    data_tmp[bar_code,]<-a$data[bar_code,]
  }
  p_tmp=Heatmap(t(as.matrix(data_tmp)), name = paste0("mat_",as.character(i+1)), 
             column_order = cell_order, 
             row_order = row_order, 
             cluster_rows = FALSE,
             cluster_columns = FALSE,show_column_names = FALSE,col = PurpleAndYellow())
  ht_list = ht_list %v% p_tmp
  i=i+1
}




ht_list = col_anno %v% p2 %v% p3 %v% p4 %v% p5 %v% p6 %v% p7 %v% p8 %v% p9

draw(ht_list,main_heatmap="mat_2")





