# Find conserved markers for the clusters. Previous name of the function was "clusterMarker_Finder"
cluster_conserved_marker_finder <- function(seuratObject, 
                                            grouping_variable = "Species", 
                                            store_dir = NULL, 
                                            store_folder = "Differentially_expressed_genes", 
                                            min_cell_for_DE = 3, 
                                            DEGtest = "bimod",
                                            lfc = 0.25, 
                                            Positive_exp = TRUE, 
                                            condition_specific_markers = FALSE,
                                            min.pct = 0.1,
                                            min.diff.pct = -Inf,
                                            verbose = TRUE,
                                            only.pos = FALSE,
                                            max.cells.per.ident = Inf,
                                            random.seed = 1,
                                            latent.vars = NULL,
                                            slot = "data",
                                            counts = numeric(),
                                            cells.1 = NULL,
                                            cells.2 = NULL,
                                            features = NULL){
  
  # merge_group_markers tells the function to find markers of each grouping variable and merge them
  # For replicates analysis we want to find markers for each replicate and find conserved markers between the replicates.
  # But in species analysis we want to find conserved markers between the species. In case if a cluster is specific to a species we want to use FindMarkers function to find markers for the species
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "Conserved_marker_grouped_by_", grouping_variable))){
    dir.create(str_c(store_dir, "/", store_folder, "/", "Conserved_marker_grouped_by_", grouping_variable), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }

  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "Conserved_marker_grouped_by_", grouping_variable, "/", "Conserved_markers_DEtest_", DEGtest))){
    dir.create(str_c(store_dir, "/", store_folder, "/", "Conserved_marker_grouped_by_", grouping_variable, "/", "Conserved_markers_DEtest_", DEGtest), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }


  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/", store_folder, "/", "Conserved_marker_grouped_by_", grouping_variable, "/", "Conserved_markers_DEtest_", DEGtest, "/")
  
  # Convert the grouping variable to a factor if it's not a factor already
  # Let's transform the metadata column into factor
  if (class(seuratObject@meta.data[[grouping_variable]]) != "factor"){
    seuratObject@meta.data[[grouping_variable]] = as.factor(seuratObject@meta.data[[grouping_variable]])
  }
  
  # Let's get the cluster idents of the seurat object
  clusterLabels = str_sort(levels(Idents(seuratObject)), numeric = TRUE)
  
  # Adding the active clustering result in the meta-data section of seurat object
  seuratObject$seurat_clusters <- Idents(seuratObject)
  
  # extract the metadata table to check how many cells are there by species in each cluster
  metaDF = seuratObject@meta.data
  
  # Let's iterate through the cluster labels
  for (i in c(1:length(clusterLabels))){
    
    # Subsetting the metadata table with the current cluster number for the active clustering solution and grouping variable
    sub_metaDF <- metaDF[metaDF[, "seurat_clusters"] == clusterLabels[i], c("seurat_clusters", grouping_variable)]
    
    # Minimum number of cells to be present in a cluster for each group (dataset/replicate/species) to consider them while looking for differentially expressed genes
    temp_groups = names(which(table(sub_metaDF[[grouping_variable]]) < min_cell_for_DE))
    
    # Creating an empty vector to store the cells to remove from the data
    cellsToRemove = c()
    
    # In this step, we go over each replicate (in case grouping variable is set to Datasets) or species (species as grouping variable)
    
    # If there are less than min_cell_for_DE (3) in more than one datasets or groups
    if (length(temp_groups) > 1){
      for (j in c(1:length(temp_groups))){
        temp_cells = rownames(sub_metaDF[sub_metaDF[, grouping_variable] == temp_groups[j], ])
        cellsToRemove = c(cellsToRemove, temp_cells)
      }
    } else {
      temp_cells = rownames(sub_metaDF[sub_metaDF[, grouping_variable] == temp_groups, ])
      cellsToRemove = c(cellsToRemove, temp_cells)
    }
    
    # get all the replicate/data/experiment ids
    temp_grouping_factors = levels(sub_metaDF[, grouping_variable])
    
    # get all the replicate/data/experiment that have more cells than the threshold (min_cell_for_DE)
    temp_remaining_factors = setdiff(temp_grouping_factors, temp_groups)
    
    # Are the replicates or datasets coming from one or more species (but this can be done only because we named the replicates in the same manner; Hirsuta (OX1, OX2, OX3)) 
    temp_remaining_groups = unique(gsub("[^A-Za-z]", "",temp_remaining_factors))
    # So, if there is only replicates of one species that have enough cell to perform conserved markers function, then we can basically perform FindMarkers function
    # In case, some of the replicates don't pass the min_cell_for_DE threshould but few replicates from both species meet the criteria, we can perform FindConservedMarkers function
    
    # Remove the cells from the total cell population to perform the differential expression analysis
    cells_to_keep = colnames(seuratObject)[!colnames(seuratObject) %in% cellsToRemove]
    
    # Creating a temporary seurat object to perform DE analysis after removing the cells
    temp_seuratObject = subset(seuratObject, cells = cells_to_keep)
    
    # The if condition checks if we have more than 1 species or not. If there are two species or conditions this condition will run the conserved markers function.
    # This condition will be implemented in case of species or conditions are there
    if (length(temp_remaining_groups) >= 2){
      # use the temp_seuratObject object to find DEGs
      #print("FindConservedMarkers function from the first if statement (for more than two grouping variable).")
      CM = FindConservedMarkers(temp_seuratObject, ident.1 = clusterLabels[i], grouping.var = grouping_variable, test.use = DEGtest, logfc.threshold = lfc, only.pos = Positive_exp)
      do.call("<-", list(str_c("Cluster_", clusterLabels[i], "_CM"), CM))
      save(list = str_c("Cluster_", clusterLabels[i], "_CM"), file = str_c(temp_dir, "Cluster_", clusterLabels[i], "_positive_conserved_marker.RData"))
    }
    
    # The else condition is applied when we have one condition or species
    # This condition will be applied when we do replicate comparison of a single species or condition.
    else {
      
      # The if condition checks  if there are two or more factors in a group. To find conserved markers between replicates or datasets this condition executed
      # condition_specific_markers tells that we are interested in finding conserved markers and not a factor or replicate specific marker
      # This condition is applied when we want to find conserved markers between replicates 
      
      if ((length(temp_remaining_factors) >= 2) & (condition_specific_markers == FALSE)){
        
        #print("FindConservedMarkers function (3rd one) from the 2nd if statement (one species, two or more replicate).")
        CM = FindConservedMarkers(temp_seuratObject, ident.1 = clusterLabels[i], grouping.var = grouping_variable, test.use = DEGtest, logfc.threshold = lfc, only.pos = Positive_exp)
        do.call("<-", list(str_c("Cluster_", clusterLabels[i], "_CM"), CM))
        save(list = str_c("Cluster_", clusterLabels[i], "_CM"), file = str_c(temp_dir, "Cluster_", clusterLabels[i], "_positive_conserved_marker.RData"))
      } 
      
      # the else if condition checks if there are two or more factors in a group, and we want to find markers for the whole group or cluster as the cluster is already species or a group specific.
      # This condition is applied when we have either a single replicate in a cluster or we want to find markers for whole cluster if the cluster is species-specific.
      # This is used if we two or more conditions, and if a cluster is species specific we can identify markers for that cluster
      else if (((length(temp_remaining_factors) >= 2) & (condition_specific_markers == TRUE)) | (length(temp_remaining_factors) == 1)){
        
        #print("FindMarkers function (2nd one) from the 2nd if statement (one species, one group).")
        CM = FindMarkers(temp_seuratObject, ident.1 = clusterLabels[i], test.use = DEGtest, logfc.threshold = lfc, only.pos = Positive_exp)
        do.call("<-", list(str_c("Cluster_", clusterLabels[i], "_M"), CM))
        save(list = str_c("Cluster_", clusterLabels[i], "_M"), file = str_c(temp_dir, "Cluster_", clusterLabels[i], "_positive_marker.RData"))
        cat(str_c("Cluster ", clusterLabels[i], " specific to ", setdiff(names(table(sub_metaDF[, grouping_variable])), temp_groups), "."), file = str_c(temp_dir, "CC_", clusterLabels[i], "_positive_markers.txt"))
      }
    }
  }
}
