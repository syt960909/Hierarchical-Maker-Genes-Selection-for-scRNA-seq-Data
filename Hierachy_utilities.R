#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Take the MAtrix behind the heatmap generated by Seurat DoHeatmap function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GetHeatmapMatrix <- function(
    object,
    features = NULL,
    cells = NULL,
    group.by = 'ident',
    group.bar = FALSE,
    group.colors = NULL,
    disp.min = -2.5,
    disp.max = NULL,
    slot = 'scale.data',
    assay = NULL,
    label = TRUE,
    size = 5.5,
    hjust = 0,
    angle = 45,
    raster = TRUE,
    draw.lines = TRUE,
    lines.width = NULL,
    group.bar.height = 0.02,
    combine = TRUE
) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- (x = unique(x = features))
  disp.max <- disp.max %||% ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  # make sure features are present
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(x = features) == 0) {
      stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot,
            " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(
    object = object,
    slot = slot)[features, cells, drop = FALSE])))
  
  return(data)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Calculate the heatmap score: Average value of the diagonal blocks - k * Average of the off-diagonal blocks
#k equals to 5 by default
#Input:
#1.object: Seurat Object
#2.data: the heatmap matrix getting from GetHeatmapMatrix() function
#3.top10: the top10 marker genes corresponding to each cluster within the current level
#4.my_levels: The lineages of current level
#Output:
#Heatmap score
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HeatmapScore <- function(
    object,
    data,
    top10,
    my_levels
    
) {
  heatmap_gene<-variable.names(data)
  HeatmapScore=0
  for (i in seq_along(my_levels)) {
    tmp_gene_diagonal_sum=0
    tmp_gene_diagonal_count=0
    tmp_gene_off_diagonal_sum=0
    tmp_gene_off_diagonal_count=0
    tmp_level_gene<-heatmap_gene[heatmap_gene %in% top10$gene[top10$cluster==my_levels[i]]]
    ###TODO:make sure the mean of diagonal and off diagonal blocks are calculated one by one
    for (j in seq_along(tmp_level_gene)){
      tmp_off_diagonal_block_max = -Inf
      tmp_block_diagonal<-data[[tmp_level_gene[j]]][rownames(data) %in% colnames(object)[object@active.ident==my_levels[i]]]
      tmp_gene_diagonal_sum = tmp_gene_diagonal_sum + sum(tmp_block_diagonal)
      tmp_gene_diagonal_count = tmp_gene_diagonal_count + length(tmp_block_diagonal)
      HeatmapScore = HeatmapScore + (tmp_gene_diagonal_sum)/(tmp_gene_diagonal_count)
      for (tmp_off_type in seq_along(my_levels[-i])){
        tmp_block_off_diagonal<-data[[tmp_level_gene[j]]][(rownames(data) %in% colnames(object)[object@active.ident==my_levels[-i][tmp_off_type]])]
        tmp_gene_off_diagonal_sum = tmp_gene_off_diagonal_sum + sum(tmp_block_off_diagonal)
        tmp_gene_off_diagonal_count = tmp_gene_off_diagonal_count + length(tmp_block_off_diagonal)
        #HeatmapScore = HeatmapScore - (tmp_gene_off_diagonal_sum)/(tmp_gene_off_diagonal_count)
        if((tmp_gene_off_diagonal_sum)/(tmp_gene_off_diagonal_count)>tmp_off_diagonal_block_max){
          tmp_off_diagonal_block_max<-(tmp_gene_off_diagonal_sum)/(tmp_gene_off_diagonal_count)
        }
      }
      HeatmapScore = HeatmapScore - 5*tmp_off_diagonal_block_max
      
    }
    #HeatmapScore = HeatmapScore + (tmp_gene_diagonal_sum)/(tmp_gene_diagonal_count)-(tmp_gene_off_diagonal_sum)/(tmp_gene_off_diagonal_count)
    
    heatmap_gene<-heatmap_gene[-which(heatmap_gene %in% top10$gene[top10$cluster==my_levels[i]])]
  }
  HeatmapScore = HeatmapScore/length(my_levels)##remove or not?
  
  return(HeatmapScore)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Find the best pairwise merge for the current level lineages
#Input:
#1.object: Seurat object
#2.best_score: heatmap score of the previous optimal merge
#Output:
#seurat object annotated by the optimal pairwise merged cluster labels
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MergeClusters <- function(object,best_score=-Inf) 
{
  ori_ident<-as.character(object@active.ident)
  ori_best_score<-best_score
  current_level<-levels(object)
  for (i in seq_along(current_level)){
    for (j in seq_along(current_level)[-(1:i)]){
      cell_type_1 = current_level[i]
      #print(cell_type_1)
      cell_type_2 = current_level[j]
      #print(cell_type_2)
      new_merge<-paste(cell_type_1,';',cell_type_2,sep = '')
      tmp_new_label<-ori_ident
      tmp_new_label[which(tmp_new_label %in% c(cell_type_1,cell_type_2))]<-new_merge
      Idents(object)<-tmp_new_label
      ### TODO:BuildClusterTree()?
      my_levels<-levels(object)
      #print(my_levels)
      object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      top10 <- object.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
      data<-GetHeatmapMatrix(object,features = top10$gene)
      
      score<-HeatmapScore(object,data,top10,my_levels)
      if (score>best_score){
        best_score<-score
        object@meta.data$best_merge_label<-tmp_new_label
        object$best_score<-best_score
        
      }
    }
  }
  if(best_score!=ori_best_score){Idents(object)<-object@meta.data$best_merge_label}
  
  
  return(object)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###write a loop to merge till length(levels)==2-->find the best merge at one node
#Input:
#object: Seurat object with optimal merge annotation for previous level
#Output:
#object: Seurat object with optimal merge annotation for current level
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FindNodeBestMerge <- function(object) 
{
  while (length(levels(object))>2) {
    best_score<-object$best_score[1]
    best_score_ori<-object$best_score[1]
    ptm <- proc.time()
    object<-MergeClusters(object,best_score = best_score)
    print('time for one merge')
    print(ptm - proc.time())
    if(object$best_score[1]>=best_score){
      best_score<-object$best_score[1]
      node_best_merge_label<-object@meta.data$best_merge_label
      print(levels(object))
    }
    if (object$best_score[1]==best_score_ori){break}###remove?
    
    Idents(object)<-node_best_merge_label
    object$best_score<-best_score
    
  }
  
  
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Write a loop to finish merging all the cell types for the child nodes
#Input:
#object:seurat object without hierarchical annotation
#Output:
##object_list:a list of seurat object with hierarchical annotation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BuildHierarchicalMap <- function(object) 
{
  object_list<-list(object)
  object_list_tmp<-list(object)
  
  while (length(object_list_tmp)!=0){
    object<-object_list_tmp[[1]]
    ptm_bm <- proc.time()
    object<-FindNodeBestMerge(object)
    print('time for best merge')
    print(ptm_bm <- proc.time())
    object_list_tmp[[1]]$best_score<-object$best_score[1]
    ####TODO:check the loop condition
    print(levels(object_list_tmp[[1]]))
    print(levels(object))
    if(length(levels(object))==length(levels(object_list_tmp[[1]]))){ 
      print('done for this node')
      object_list_tmp<-object_list_tmp[-1]
      for (tmp in object_list_tmp){
        print(levels(tmp))
      }
      
      next
    }
    object_list<-c(object_list,object)
    object_levels<-levels(object)
    for (i in seq_along(object_levels)) 
    {
      sub_cell_type_group<-strsplit(object_levels[i], ";")[[1]]
      if(length(sub_cell_type_group)<=2){
        for (tmp in object_list_tmp){
          print(levels(tmp))
        }
        
        next
      }
      ###Now we have a list of seurat object that each element of it may need further mergeing(each element can be seen as a child node)
      object_sub<-subset(x = object_list[[1]], idents = strsplit(object_levels[i], ";")[[1]])
      ###TODO:a better way to define the initial best score of the subgroup?
      my_levels<-levels(object_sub)
      #print(my_levels)
      object.markers <- FindAllMarkers(object_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      top10 <- object.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
      data_sub<-GetHeatmapMatrix(object_sub,features = top10$gene)
      
      score<-HeatmapScore(object_sub,data_sub,top10,my_levels)
      object_sub$best_score<-score
      object_list_tmp<-c(object_list_tmp,object_sub)
      
    }
    object_list_tmp<-object_list_tmp[-1]
    print('test 1')
    for (tmp in object_list_tmp)
    {
      print(levels(tmp))
    }
  }
  return(object_list)
}


PurpleAndYellow <- function(k = 50) {
  return(CustomPalette(low = "magenta", high = "yellow", mid = "black", k = k))
}


DoHeatmap_matrix <- function(
    object,
    features = NULL,
    cells = NULL,
    group.by = 'ident',
    group.bar = TRUE,
    group.colors = NULL,
    disp.min = -2.5,
    disp.max = NULL,
    slot = 'scale.data',
    assay = NULL,
    label = TRUE,
    size = 5.5,
    hjust = 0,
    angle = 45,
    raster = TRUE,
    draw.lines = FALSE,
    lines.width = NULL,
    group.bar.height = 0.02,
    combine = TRUE
) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  # make sure features are present
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(x = features) == 0) {
      stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot,
            " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(
    object = object,
    slot = slot)[features, cells, drop = FALSE])))
  data_ = data
  object <- suppressMessages(expr = StashIdent(object = object, save.name = 'ident'))
  group.by <- group.by %||% 'ident'
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  # group.use <- switch(
  #   EXPR = group.by,
  #   'ident' = Idents(object = object),
  #   object[[group.by, drop = TRUE]]
  # )
  # group.use <- factor(x = group.use[cells])
  plots <- vector(mode = 'list', length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    #print('########################group.use####################')
    #print(group.use)
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    #print('######################cells############################')
    #print(cells)
    names(x = group.use) <- cells
    if (draw.lines) {
      # create fake cells to serve as the white lines, fill with NAs
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 0.0025)
      placeholder.cells <- sapply(
        X = 1:(length(x = levels(x = group.use)) * lines.width),
        FUN = function(x) {
          return(RandomName(length = 20))
        }
      )
      placeholder.groups <- rep(x = levels(x = group.use), times = lines.width)
      group.levels <- levels(x = group.use)
      names(x = placeholder.groups) <- placeholder.cells
      group.use <- as.vector(x = group.use)
      names(x = group.use) <- cells
      group.use <- factor(x = c(group.use, placeholder.groups), levels = group.levels)
      na.data.group <- matrix(
        data = NA,
        nrow = length(x = placeholder.cells),
        ncol = ncol(x = data.group),
        dimnames = list(placeholder.cells, colnames(x = data.group))
      )
      data.group <- rbind(data.group, na.data.group)
    }
    lgroup <- length(levels(group.use))
    # plot <- SingleRasterMap(
    #   data = data.group,
    #   raster = raster,
    #   disp.min = disp.min,
    #   disp.max = disp.max,
    #   feature.order = features,
    #   cell.order = names(x = sort(x = group.use)),
    #   group.by = group.use
    # )
    # plot <- plot + theme(line = element_blank())
    # plots[[i]] <- plot
    my_list <- list("data" = data.group, "disp.min" = disp.min, 'disp.max' = disp.max,
                    'feature.order' = features,
                    'cell.order' = names(x = sort(x = group.use)),
                    'group.by' = group.use
    )
  }
  
  return(my_list)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Generate hierarchical tree based on BuildClusterTree function in Seurat, which constructs a phylogenetic tree relating the 'average' cell from each identity class.
#Tree is estimated based on a distance matrix constructed in either gene expression space or PCA space(default).

#Input:
#object:seurat object without hierarchical annotation
#Output:
##a list with components:
#1.merge:an n-1 by 2 matrix. Row i of merge describes the merging of clusters at step i of the clustering. If an element j in the row is negative, then observation -j was merged at this stage. If j is positive then the merge was with the cluster formed at the (earlier) stage j of the algorithm. Thus negative entries in merge indicate agglomerations of singletons, and positive entries indicate agglomerations of non-singletons.
#2.height:a set of n-1 real values (non-decreasing for ultrametric trees). The clustering height: that is, the value of the criterion associated with the clustering method for the particular agglomeration.
#3.order:a vector giving the permutation of the original observations suitable for plotting, in the sense that a cluster plot using this ordering and matrix merge will not have crossings of the branches.
#4.labels:labels for each of the objects being clustered.
#5.call:the call which produced the result.
#6.method:the cluster method that has been used.
#7.dist.method:the distance that has been used to create d (only returned if the distance object has a “method” attribute).
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BuildClusterTree_predefine <- function(
    object,
    assay = NULL,
    features = NULL,
    dims = NULL,
    reduction = "pca",
    graph = NULL,
    slot = 'data',
    reorder = FALSE,
    reorder.numeric = FALSE,
    verbose = TRUE
) {
  if (!PackageCheck('ape', error = FALSE)) {
    stop(cluster.ape, call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = graph)) {
    idents <- levels(x = object)
    nclusters <- length(x = idents)
    data.dist <- matrix(
      data = numeric(length = 1L),
      nrow = nclusters,
      ncol = nclusters,
      dimnames = list(idents, idents)
    )
    graph <- object[[graph]]
    cxi <- CellsByIdentities(object = object)
    cpairs <- na.omit(object = unique(x = t(x = apply(
      X = expand.grid(1:nclusters, 1:nclusters)[, c(2, 1)],
      MARGIN = 1,
      FUN = function(x) {
        if (length(x = x) == length(x = unique(x = x))) {
          return(sort(x = x))
        }
        return(c(NA, NA))
      }
    ))))
    if (verbose) {
      pb <- txtProgressBar(style = 3, file = stderr())
    }
    for (i in 1:nrow(x = cpairs)) {
      i1 <- cpairs[i, ][1]
      i2 <- cpairs[i, ][2]
      graph.sub <- graph[cxi[[idents[i1]]], cxi[[idents[i2]]]]
      d <- mean(x = graph.sub)
      if (is.na(x = d)) {
        d <- 0
      }
      data.dist[i1, i2] <- d
      if (verbose) {
        setTxtProgressBar(pb = pb, value = i / nrow(x = cpairs))
      }
    }
    if (verbose) {
      close(con = pb)
    }
    diag(x = data.dist) <- 1
    data.dist <- dist(x = data.dist)
  } else if (!is.null(x = dims)) {
    my.lapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
    embeddings <- Embeddings(object = object, reduction = reduction)[, dims]
    data.dims <- my.lapply(
      X = levels(x = object),
      FUN = function(x) {
        cells <- WhichCells(object = object, idents = x)
        if (length(x = cells) == 1) {
          cells <- c(cells, cells)
        }
        temp <- colMeans(x = embeddings[cells, ])
      }
    )
    data.dims <- do.call(what = 'cbind', args = data.dims)
    colnames(x = data.dims) <- levels(x = object)
    data.dist <- dist(x = t(x = data.dims))
  } else {
    features <- features %||% VariableFeatures(object = object)
    features <- intersect(x = features, y = rownames(x = object))
    data.avg <- AverageExpression(
      object = object,
      assays = assay,
      features = features,
      slot = slot,
      verbose = verbose
    )[[1]]
    data.dist <- dist(x = t(x = data.avg[features, ]))
  }
  #data.tree <- ape::as.phylo(x = hclust(d = data.dist))
  data.tree <- hclust(d = data.dist)
  return(data.tree)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Merge lineages according to the hierarchical tree getting from BuildClusterTree_predefine() function
#Input:
#1. object: Seurat object without hierarchical annotation
#2. predefined_merge_list:the hierarchy predifined by BuildClusterTree_predefine
#3. best_score: the heatmap score of the seurat object without hierarchical annotation
#Output:
#object: Seurat object with the optimal grouped clusters annotation in terms of the clusters within it.

MergeClusters_predefine <- function(object,predefined_merge_list,best_score=-Inf) 
{
  ori_ident<-as.character(object@active.ident)
  ori_best_score<-best_score
  
  node_best_merge<-ori_ident
  for (current_level in predefined_merge_list) {
    tmp_new_label<-ori_ident
    ##############TODO###################
    for(current_merge in current_level){
      ########TODO###############
      new_merge<-paste(current_merge, collapse = ';')
      tmp_new_label[which(tmp_new_label %in% as.character(current_merge))]<-new_merge
    }
    Idents(object)<-tmp_new_label
    ### TODO:BuildClusterTree()?
    my_levels<-levels(object)
    if(length(my_levels)<2){break}
    object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    top10 <- object.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
    data<-GetHeatmapMatrix(object,features = top10$gene)
    score<-HeatmapScore(object,data,top10,my_levels)
    print(score)
    print(best_score)
    if (score<=best_score){
      Idents(object)<-node_best_merge
      print('should return and stop')
      break
      #stop
    }###remove?
    if (score>best_score){
      best_score<-score
      object@meta.data$best_merge_label<-tmp_new_label
      object$best_score<-best_score
      node_best_merge<-tmp_new_label
      
    }
  }
  
  #if(best_score!=ori_best_score){Idents(object)<-object@meta.data$best_merge_label}
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Here we provide a faster implementation based on Seurat BuildClusterTree function, which constructs a phylogenetic tree relating the 'average' cell from each identity class.
#Tree is estimated based on a distance matrix constructed in either gene expression space or PCA space(default).
#Through the BuildClusterTree, we can get a vector called "height", which is a set of n-1 real values (non-decreasing for ultrametric trees). 
#The clustering height: that is, the value of the criterion associated with the clustering method for the particular agglomeration.
#From bottom to top, we set each height value as threshold to cut the tree into groups so that we get the hierarchical lineages
#We calculate the heatmap score corresponding to lineages for each level, and stop the process when the score is not getting higher
#We keep the lineages under the stop height as the final hierarchical tree. 
#Input:
#object:seurat object without hierarchical annotation
#Output:
##object_list:a list of seurat object with hierarchical annotation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BuildHierarchicalMap_predefine <- function(object) 
{
  object_list<-list(object)
  object_list_tmp<-list(object)
  while (length(object_list_tmp)!=0){
    object<-object_list_tmp[[1]]
    #object_datatree <- BuildClusterTree_(object = object,slot = 'scale.data')
    object_datatree <- BuildClusterTree_predefine(object = object,slot = 'data')
    predefined_merge_list<-list()
    #if (length(object_datatree$height)>10){
    # tree_list<-object_datatree$height[round(length(object_datatree$height)/2):length(object_datatree$height)]
    #}
    #else{tree_list<-object_datatree$height}
    #for (i in tree_list){
    for (i in object_datatree$height){
      object_datatree_<-cutree(object_datatree,h=i)
      current_level<-list()
      for (j in names(table(object_datatree_))){
        if (table(object_datatree_)[[j]]>1){
          current_level<-c(current_level,list(names(object_datatree_[object_datatree_ %in% j])))
          print(current_level) 
          
        }
      }
      predefined_merge_list<-append(predefined_merge_list,list(current_level))
    }
    predefined_merge_list<-predefined_merge_list[1:length(predefined_merge_list)-1]
    #object<-FindNodeBestMerge_predefine(object)
    #####the code written should just be a substitute of FindNodeBestMerge()
    
    object<-MergeClusters_predefine(object,predefined_merge_list,best_score=object$best_score[[1]])
    
    print('###########next step#############')
    object_list_tmp[[1]]$best_score<-object$best_score[1]
    ####TODO:check the loop condition
    print(levels(object))
    if(length(levels(object))==length(levels(object_list_tmp[[1]]))){ 
      print('done for this node')
      object_list_tmp<-object_list_tmp[-1]
      next
    }
    
    object_list<-c(object_list,object)
    object_levels<-levels(object)
    for (i in seq_along(object_levels)) 
    {
      sub_cell_type_group<-strsplit(object_levels[i], ";")[[1]]
      if(length(sub_cell_type_group)<=2){
        next
      }
      ###Now we have a list of seurat object that each element of it may need further merging(each element can be seen as a child node)
      object_sub<-subset(x = object_list[[1]], idents = strsplit(object_levels[i], ";")[[1]])
      ###TODO:a better way to define the initial best score of the subgroup?
      my_levels<-levels(object_sub)
      print('#########sub_level##########')
      print(my_levels)
      object.markers <- FindAllMarkers(object_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      top10 <- object.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
      data_sub<-GetHeatmapMatrix(object_sub,features = top10$gene)
      
      score<-HeatmapScore(object_sub,data_sub,top10,my_levels)
      object_sub$best_score<-score
      object_list_tmp<-c(object_list_tmp,object_sub)
      
    }
    object_list_tmp<-object_list_tmp[-1]
    for (tmp in object_list_tmp)
    {
      print('#######################')
      print(levels(tmp))
    }
  }
  return(object_list)
}
