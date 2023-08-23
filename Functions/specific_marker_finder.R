specific_marker_finder <- function(DEG_file_dir, 
                                   pct_detection = 0.1, 
                                   pct_diff = 0.3,
                                   store_dir, 
                                   store_folder = "Cluster_specific_markers",
                                   store.marker.file = FALSE,
                                   return_markers_list = FALSE,
                                   marker_file_name = NULL){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # Initializing an empty list to store all the DEG files
  temp_deg_list = list()
  
  # Check the files and load them
  if (class(DEG_file_dir) == "character"){
    
    DEG_files = str_sort(list.files(DEG_file_dir, pattern = "Cluster_"), numeric = TRUE)
    
    for (i in c(1:length(DEG_files))){
      
      SM = loadRData(str_c(DEG_file_dir, "/", DEG_files[i]))
      
      temp_deg_list[[i]] <- SM
      
      names(temp_deg_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]), "_SM")
    }
  }
  
  else if (class(DEG_file_dir) != "character"){
    
    if (class(DEG_file_dir) != "list"){
      temp_deg_list = list(DEG_file_dir)
      
      if (length(temp_deg_list) == 1){
        
        names(temp_deg_list) <- "Cluster_specific_markers"
      }
    }
  }
  
  # Let's create a list to store the specific marker files to the list
  temp_list = list()
  
  # Get the names of the list items or DEG files
  temp_names = names(temp_deg_list)
  
  
  for (i in c(1:length(temp_deg_list))){
    
    # load the DEG file
    Cluster_DEGs = temp_deg_list[[temp_names[i]]]
    
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
    
    
    # checks if the user want to save the markers file in the directory or not
    if (store.marker.file == TRUE){
      
      if (grepl("[[:digit:]]", temp_names[i]) != TRUE) {
        
        if (missing(marker_file_name)) {
          
          write.csv(Cluster_specific_DEGs, file = "Cluster_specific_markers.csv", row.names = TRUE)
        }
        
        else {
          
          write.csv(Cluster_specific_DEGs, file = str_c(marker_file_name, ".csv"), row.names = TRUE)
          
        }
      }
      
      else {
        # Save the file
        write.csv(Cluster_specific_DEGs, file = str_c(temp_dir, "Cluster_", parse_number(temp_names[i]), "_specific_markers.csv"), row.names = TRUE)
      }
    }
    
    else {
      Cluster_specific_DEGs
    }
    
    
    # Storing the specific marker file to the list
    temp_list[[i]] <- Cluster_specific_DEGs
    
    if (grepl("[[:digit:]]", temp_names[i]) != TRUE) {
      
      names(temp_list)[i] <- str_c("Cluster_specific_markers")
      
    }
    
    else {
      names(temp_list)[i] <- str_c("Cluster_", parse_number(temp_names[i]))
    }
  }
  
  if (return_markers_list == TRUE){
    return(temp_list)
  }
}
