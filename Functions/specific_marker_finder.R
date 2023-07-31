specific_marker_finder <- function(DEG_file_dir, pct_detection = 0.1, pct_diff = 0.3){
  
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
    
    
    temp_list[[i]] <- Cluster_specific_DEGs
    names(temp_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
  }
  
    return(temp_list)

}
