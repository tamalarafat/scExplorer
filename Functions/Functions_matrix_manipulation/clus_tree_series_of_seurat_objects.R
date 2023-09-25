# For a series of seurat objects
clus_tree_series_of_seurat_objects <- function(seurat_object_dir, 
                                               file_name_pattern, 
                                               marker_file, 
                                               query_pattern, 
                                               gene_ID_column, 
                                               gene_name_column = NULL, 
                                               store_dir = NULL, 
                                               store_folder = "Tree_of_cell_clusters_graph_based_clustering"){
  
  # Lets get the saved file names
  seuratObjects = str_sort(list.files(path = seurat_object_dir, pattern = file_name_pattern), numeric = TRUE)
  
  for (i in c(1:length(seuratObjects))){
    # Creating the efile name with path location
    file_name = str_c(seurat_object_dir, "/", seuratObjects[i])
    
    # Loading the liger object and storing it to a temporary variable
    integrated.data = loadRData(file_name) # Load the RData file and assign it to a variable using the function loaadRData
    
    DefaultAssay(integrated.data) <- "RNA"
    
    clus_tree_of_cell_clusters(seuratObject = integrated.data, marker_file = marker_file, query_pattern = query_pattern, gene_ID_column = gene_ID_column, gene_name_column = gene_name_column, store_dir)
  }
}
