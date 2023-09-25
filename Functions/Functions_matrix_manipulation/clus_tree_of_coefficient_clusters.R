# Generate cluster tree of the quantile normalized cell clusters
clus_tree_of_coefficient_clusters <- function(seuratObject, 
                                              query_pattern, 
                                              graph_based_clustering = FALSE, 
                                              store_dir = NULL, 
                                              store_folder = NULL,
                                              figure_name = "Arrangement_of_cell_clusters.png"){
  
  # Let's create a directory to store the GO annotation results for different factorization
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  # Let's create a directory to store the GO annotation results
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  if (graph_based_clustering == TRUE){
    
    # Generate the cluster tree for the cell clustering for different resolution parameter
    p <- clustree(seuratObject, prefix = query_pattern, node_text_size = 8, node_size = 14) + 
      labs(colour = "Resolution\nvalue", edge_colour = "Count", edge_alpha = "Proportion") + 
      theme(legend.key.size = unit(2,"line"),
            legend.text = element_text(size = 20, face = "bold"),
            legend.title = element_text(size = 22, colour = "black", face = "bold"))
    ggsave(filename = str_c(temp_dir, figure_name), plot = p, width = 48, height = 22, dpi = 300)
    
  }
  
  else {
    mat = seuratObject@meta.data
    mat = mat[ , str_sort(colnames(mat)[str_detect(colnames(mat), query_pattern)], numeric = TRUE)]
    
    # dividing the columns to specific length
    temp = seq(1, ncol(mat), 10)
    
    tempNames = str_sort(colnames(mat), numeric = TRUE)
    
    for (i in c(1:(length(temp)-1))){
      tempMat = mat[ ,c(temp[i]:(temp[(i+1)] - 1))]
      
      p <- clustree(tempMat, prefix = unique(sub("(\\d)+", "", tempNames)), node_text_size = 8, node_size = 14) + 
        labs(colour = "Factorization\nK", edge_colour = "Count", edge_alpha = "Proportion") + 
        theme(legend.key.size = unit(2,"line"),
              legend.text = element_text(size = 20, face = "bold"),
              legend.title = element_text(size = 22, colour = "black", face = "bold"))
      
      ggsave(filename = str_c(temp_dir, "Part_", i, "_Tree_of_cell_clusters.png"), plot = p, width = 48, height = 22, dpi = 300)
    }
    
    temp_df = mat[ , temp[length(temp)]:ncol(mat)]
    p <- clustree(temp_df, prefix = unique(sub("(\\d)+", "", tempNames)), node_text_size = 8, node_size = 14) + 
      labs(colour = "Factorization\nK", edge_colour = "Count", edge_alpha = "Proportion") + 
      theme(legend.key.size = unit(2,"line"),
            legend.text = element_text(size = 20, face = "bold"),
            legend.title = element_text(size = 22, colour = "black", face = "bold"))
    
    ggsave(filename = str_c(temp_dir, "Part_", i+1, "_Tree_of_cell_clusters.png"), plot = p, width = 48, height = 22, dpi = 300)
    
    # Let's create a basis matrix with only a few components
    tempVec = as.numeric(str_extract(tempNames, "(\\d)+"))
    miniBatch = unique(c(1, which(tempVec %in% tempVec[tempVec%%5 == 0])))
    
    if (max(miniBatch) > ncol(mat)){
      miniBatch[which.max(miniBatch)] = ncol(mat)
    }
    
    miniCoef = mat[, miniBatch]
    
    p <- clustree(miniCoef, prefix = unique(sub("(\\d)+", "", tempNames)), node_text_size = 8, node_size = 14) + 
      labs(colour = "Factorization\nK", edge_colour = "Count", edge_alpha = "Proportion") + 
      theme(legend.key.size = unit(2,"line"),
            legend.text = element_text(size = 20, face = "bold"),
            legend.title = element_text(size = 22, colour = "black", face = "bold"))
    
    ggsave(filename = str_c(temp_dir, "Tree_of_cell_clusters.png"), plot = p, width = 48, height = 22, dpi = 300)
  }
}

