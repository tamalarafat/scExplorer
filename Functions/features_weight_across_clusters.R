features_weight_across_clusters <- function(seuratObject, 
                                            features_list, 
                                            search_max, 
                                            include_cell_count, 
                                            store_dir = NULL, 
                                            store_folder = "Duplicated_genes_details", 
                                            file_name = "Duplicated_genes_distribution",
                                            return_output = FALSE) {
  
  # Creating necessary storing space to store the results
  
  # Path to store the output
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  # Directory to store the output
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # The path of the output folder
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  
  # The input features set can be a vector or data.frame 
  if (features_list == "character") {
    genes_file = data.frame(gene_ID = features_list)
  }
  
  else if (features_list == "data.frame") {
    genes_file = features_list
  }
  
  
  # Checking if all the genes of the provided features set present in the Seurat object. 
  # We will keep only those genes that are present, and will only search on the Seurat object with these genes.
  df_present = genes_file[genes_file$gene_ID %in% intersect(genes_file$gene_ID, rownames(seuratObject)), ]
  
  # Setting rownames of the dataframe with the gene IDs
  rownames(df_present) = df_present$gene_ID
  
  # Adding columns to the data.frame to store the desired information - cluster ID with max number of expressing cells, number of cells expressing the gene in the clusters, etc..
  df_present[, c(str_c("cluster_ID_max_", english(c(1:search_max))), str_c("cell_count_max_", english(c(1:search_max))))] = ""
  
  # Genes that are absent in the Seurat object
  df_absent = genes_file[!genes_file$gene_ID %in% intersect(genes_file$gene_ID, rownames(seuratObject)), ]
  
  # Adding columns to the data.frame to store the desired information - "Absent"
  df_absent[, c(str_c("cluster_ID_max_", english(c(1:search_max))), str_c("cell_count_max_", english(c(1:search_max))))] = "NA"
  
  # Setting rownames of the dataframe with the gene IDs
  rownames(df_absent) = df_absent$gene_ID
  
  # Get the active cluster labels in the Seurat object
  cluster_labels = str_sort(levels(Idents(seuratObject)), numeric = TRUE)
  
  for (geneID in c(1:nrow(df_present))){
    
    # Lets create a empty vector to store the count information
    temp = list()
    
    # iterate over the cluster labels to get the number of cells in each cluster expressing the given gene
    for (i in c(1:length(cluster_labels))){
      
      # get the cell IDs belonging to each cluster
      cell_names = WhichCells(seuratObject, idents = cluster_labels[i])
      
      # Number of cells expressing the gene for each cluster
      cell_count = sum(GetAssayData(seuratObject, assay = "RNA", slot = "counts")[df_present$gene_ID[geneID], cell_names] != 0)
      
      # For each cluster, adding the information to the list 
      temp[[i]] = data.frame(cluster_ID = cluster_labels[i], expressing_cells = cell_count)
    }
    
    # data.frame containing the information generated from the iteration
    df_count_holder = do.call(rbind.data.frame, temp)
    
    # The search max specifies how many cluster IDs with max counts to keep
    df_count_holder = df_count_holder[order(df_count_holder$expressing_cells, decreasing = TRUE)[1:search_max], ]
    
    # Adding the information to the data.frame containing the gene IDs
    df_present[df_present$gene_ID[geneID], c(str_c("cluster_ID_max_", english(c(1:search_max))), str_c("cell_count_max_", english(c(1:search_max))))] = c(df_count_holder$cluster_ID[c(1:search_max)], df_count_holder$expressing_cells[c(1:search_max)])
    
  }
  
  genes_file = rbind.data.frame(df_present, df_absent)
  
  if (include_cell_count == TRUE) {
    
    write.csv(genes_file, file = str_c(temp_dir, file_name, ".csv"), row.names = FALSE)
  } else {
    
    write.csv(genes_file[, -grep(pattern = "cell_count_max", x = colnames(genes_file))], file = str_c(temp_dir, file_name, ".csv"), row.names = FALSE)
  }
  
  
  if (return_output == TRUE) {return(genes_file)}
}

