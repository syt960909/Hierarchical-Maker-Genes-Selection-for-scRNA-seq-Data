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
      top10 <- object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
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




###write a loop to merge till length(levels)==2-->find the best merge at one node

FindNodeBestMerge <- function(object) 
{
  while (length(levels(object))>2) {
    best_score<-object$best_score[1]
    best_score_ori<-object$best_score[1]
    object<-MergeClusters(object,best_score = best_score)
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



### Write a loop to finish merging all the cell types for the child nodes
BuildHierarchicalMap <- function(object) 
{
  object_list<-list(object)
  object_list_tmp<-list(object)
  
  while (length(object_list_tmp)!=0){
    object<-object_list_tmp[[1]]
    print('========test 4=====================')
    object<-FindNodeBestMerge(object)
    print('========test 5=====================')
    object_list_tmp[[1]]$best_score<-object$best_score[1]
    ####TODO:check the loop condition
    print('test 2')
    print(levels(object_list_tmp[[1]]))
    print(levels(object))
    if(length(levels(object))==length(levels(object_list_tmp[[1]]))){ 
      print('done for this node')
      object_list_tmp<-object_list_tmp[-1]
      print('========test 7=====================')
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
      print('========test 8=====================')  
      ###Now we have a list of seurat object that each element of it may need further mergeing(each element can be seen as a child node)
      object_sub<-subset(x = object_list[[1]], idents = strsplit(object_levels[i], ";")[[1]])
      ###TODO:a better way to define the initial best score of the subgroup?
      my_levels<-levels(object_sub)
      #print(my_levels)
      object.markers <- FindAllMarkers(object_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      top10 <- object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
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
