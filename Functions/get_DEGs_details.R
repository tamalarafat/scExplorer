get_DEGs_details <- function(DEG_file_dir, gene_description_file, store_dir, store_folder = "Cluster_DEGs_and_Specific_markers", return_markers_list = TRUE){
 
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  DEG_files = str_sort(list.files(DEG_file_dir, pattern = "Cluster_"), numeric = TRUE)
  
  CH_temp_cluster_DEG_list = list()
  CH_specific_cluster_marker_list = list()
  AT_temp_cluster_DEG_list = list()
  AT_specific_cluster_marker_list = list()
  
  for (i in c(1:length(DEG_files))){
    
    # Let's create a folder to store all the markers file
    if (!dir.exists(str_c(temp_dir, "Cluster_", parse_number(DEG_files[i])))){
      dir.create(str_c(temp_dir, "Cluster_", parse_number(DEG_files[i])), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    temp_store_dir = str_c(temp_dir, "Cluster_", parse_number(DEG_files[i]), "/")
    
    # load the DEG file
    Cluster_DEGs = loadRData(str_c(DEG_file_dir, "/", DEG_files[i]))
    
    # Prepare the file of cluster differentially expressed gene
    Cluster_DEGs = Cluster_DEGs[order(Cluster_DEGs$pct.2, decreasing = FALSE), ]
    Cluster_DEGs$gene_ID = rownames(Cluster_DEGs)
    
    # Save the file
    write.csv(Cluster_DEGs, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs.csv"), row.names = FALSE)
    
    # Write a text file containing Cardamine gene IDs
    write_lines(rownames(Cluster_DEGs), file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs_CH_ID.txt"))
    
    # Create ortho genes list and gene description file
    description = gene_description_file[rownames(gene_description_file) %in% rownames(Cluster_DEGs), ]
    
    
    # Description - All DEGs identified for the cluster
    All_DEGs_description = merge(Cluster_DEGs, description, by.x = "gene_ID", by.y = "CH_ID", all.x = TRUE)[-c(8, 9)]
    All_DEGs_description = All_DEGs_description[order(All_DEGs_description$pct.2, decreasing = FALSE), ]
    rownames(All_DEGs_description) = All_DEGs_description$gene_ID
    
    # Save the file
    write.csv(All_DEGs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs_description.csv"), row.names = FALSE)
    
    # Write the orthologues ID file
    write_lines(All_DEGs_description$AT_ID[!is.na(All_DEGs_description$AT_ID)], file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs_orthoID.txt"))
    
    # Let's get all TFs information
    All_TFs_description = All_DEGs_description[!is.na(All_DEGs_description$TF), ]
    
    # Save the file
    write.csv(All_TFs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_TF_description.csv"), row.names = FALSE)
    
    # Write the TF hirsuta ID file
    write_lines(All_TFs_description$gene_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_TF_CH_ID.txt"))
    
    # Write the TF ortho ID file
    write_lines(All_TFs_description$AT_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_TF_orthoID.txt"))
    
    # Here create a function for specific marker selection
    Cluster_specific_DEGs = Cluster_DEGs[(Cluster_DEGs$pct.2 <= 0.1) | ((Cluster_DEGs$pct.1 - Cluster_DEGs$pct.2) >= 0.3), ]
    Cluster_specific_DEGs = Cluster_specific_DEGs[order(Cluster_specific_DEGs$pct.2, decreasing = FALSE), ]
    Cluster_specific_DEGs$gene_ID = rownames(Cluster_specific_DEGs)
    
    # Save the file
    write.csv(Cluster_specific_DEGs, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers.csv"), row.names = FALSE)
    
    # Write a text file containing Cardamine gene IDs
    write_lines(rownames(Cluster_specific_DEGs), file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_CH_ID.txt"))
    
    CS_DEGs_description = merge(Cluster_specific_DEGs, description, by.x = "gene_ID", by.y = "CH_ID", all.x = TRUE)[-c(8, 9)]
    CS_DEGs_description = CS_DEGs_description[order(CS_DEGs_description$pct.2, decreasing = FALSE), ]
    rownames(CS_DEGs_description) = CS_DEGs_description$gene_ID
    
    # Save the file
    write.csv(CS_DEGs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_description.csv"), row.names = FALSE)
    
    # Write the orthologues ID file
    write_lines(CS_DEGs_description$AT_ID[!is.na(CS_DEGs_description$AT_ID)], file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_orthoID.txt"))
    
    # Description - Cluster specific TFs 
    CS_TFs_description = CS_DEGs_description[!is.na(CS_DEGs_description$TF), ]
    
    # Save the file
    write.csv(CS_TFs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_description.csv"), row.names = FALSE)
    
    # Write the TF hirsuta ID file
    write_lines(CS_TFs_description$gene_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_CH_ID.txt"))
    
    # Write the TF ortho ID file
    write_lines(CS_TFs_description$AT_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_orthoID.txt"))
    
    CH_temp_cluster_DEG_list[[i]] <- rownames(Cluster_DEGs)
    names(CH_temp_cluster_DEG_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
    
    CH_specific_cluster_marker_list[[i]] <- rownames(Cluster_specific_DEGs)
    names(CH_specific_cluster_marker_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
    
    AT_temp_cluster_DEG_list[[i]] <- All_DEGs_description$AT_ID[!is.na(All_DEGs_description$AT_ID)]
    names(AT_temp_cluster_DEG_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
    
    AT_specific_cluster_marker_list[[i]] <- CS_DEGs_description$AT_ID[!is.na(CS_DEGs_description$AT_ID)]
    names(AT_specific_cluster_marker_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
  }
  
  if (return_markers_list){
    markerList = list(cluster_all_degs_CH = CH_temp_cluster_DEG_list,
                      cluster_all_degs_ortho = AT_temp_cluster_DEG_list,
                      cluster_specific_marker_CH = CH_specific_cluster_marker_list,
                      cluster_specific_marker_ortho = AT_specific_cluster_marker_list)
    
    return(markerList)
  }
  
}
