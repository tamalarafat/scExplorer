specific_marker_finder <- function(DEG_file_dir, 
                                   pct_detection = 0.1, 
                                   pct_diff = 0.3,
                                   store_dir, 
                                   store_folder = "Cluster_specific_markers",
                                   return_markers_list = FALSE){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  DEG_files = str_sort(list.files(DEG_file_dir, pattern = "Cluster_"), numeric = TRUE)
  
  temp_list = list()
  
  for (i in c(1:length(DEG_files))){
    
    # load the DEG file
    Cluster_DEGs = loadRData(str_c(DEG_dir, "/", DEG_files[i]))
    
    # Prepare the file of cluster differentially expressed gene
    Cluster_DEGs = Cluster_DEGs[order(Cluster_DEGs$pct.2, decreasing = FALSE), ]
    Cluster_DEGs$gene_ID = rownames(Cluster_DEGs)
    
    # Let's create the selection criteria
    criteria_1 = Cluster_DEGs$pct.2 <= pct_detection
    
    criteria_2 = (Cluster_DEGs$pct.1 - Cluster_DEGs$pct.2) >= pct_diff
    # Here create a function for specific marker selection
    Cluster_specific_DEGs = Cluster_DEGs[criteria_1 | criteria_2, ]
    Cluster_specific_DEGs = Cluster_specific_DEGs[order(Cluster_specific_DEGs$pct.2, decreasing = FALSE), ]
    Cluster_specific_DEGs$gene_ID = rownames(Cluster_specific_DEGs)
    
    # Save the file
    write.csv(Cluster_specific_DEGs, file = str_c(temp_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers.csv"), row.names = FALSE)
    
    temp_list[[i]] <- Cluster_specific_DEGs
    names(temp_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
  }
  
  if (return_markers_list == TRUE){
    return(temp_list)
  }
}
