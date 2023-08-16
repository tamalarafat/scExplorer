get_conserved_DEGs_details_thaliana <- function(DEG_file_dir, 
                                                gene_description_file, 
                                                store_dir, 
                                                store_folder = "Cluster_DEGs_and_Specific_markers", 
                                                return_markers_list = TRUE){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  DEG_files = str_sort(list.files(DEG_file_dir, pattern = "Cluster_"), numeric = TRUE)
  
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
    
    # Extract all the columns in the DEG file containing pct.2 information
    temp_all_pct2 = Cluster_DEGs[ , str_detect(colnames(Cluster_DEGs), pattern = "pct.2"), ]
    
    if (class(temp_all_pct2) != "data.frame"){
      # Prepare the file of cluster differentially expressed gene
      Cluster_DEGs = Cluster_DEGs[order(Cluster_DEGs$pct.2, decreasing = FALSE), ]
      Cluster_DEGs$gene_ID = rownames(Cluster_DEGs)
      
      # Save the file
      write.csv(Cluster_DEGs, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs.csv"), row.names = FALSE)
      
      # Write a text file containing Cardamine gene IDs
      write_lines(rownames(Cluster_DEGs), file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs_AT_ID.txt"))
      
      # Create ortho genes list and gene description file
      description = gene_description_file[rownames(gene_description_file) %in% rownames(Cluster_DEGs), ]
      
      
      # Description - All DEGs identified for the cluster
      All_DEGs_description = merge(Cluster_DEGs, description, by.x = "gene_ID", by.y = "AT_ID", all.x = TRUE)[-c(8, 9)]
      All_DEGs_description = All_DEGs_description[order(All_DEGs_description$pct.2, decreasing = FALSE), ]
      rownames(All_DEGs_description) = All_DEGs_description$gene_ID
      
      # Save the file
      write.csv(All_DEGs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs_description.csv"), row.names = FALSE)
      
      # Let's get all TFs information
      All_TFs_description = All_DEGs_description[!is.na(All_DEGs_description$TF), ]
      
      # Save the file
      write.csv(All_TFs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_TF_description.csv"), row.names = FALSE)
      
      # Write the TF hirsuta ID file
      write_lines(All_TFs_description$gene_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_TF_AT_ID.txt"))
      
      # Here create a function for specific marker selection
      Cluster_specific_DEGs = Cluster_DEGs[(Cluster_DEGs$pct.2 <= 0.1) | ((Cluster_DEGs$pct.1 - Cluster_DEGs$pct.2) >= 0.3), ]
      Cluster_specific_DEGs = Cluster_specific_DEGs[order(Cluster_specific_DEGs$pct.2, decreasing = FALSE), ]
      Cluster_specific_DEGs$gene_ID = rownames(Cluster_specific_DEGs)
      
      # Save the file
      write.csv(Cluster_specific_DEGs, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers.csv"), row.names = FALSE)
      
      # Write a text file containing Cardamine gene IDs
      write_lines(rownames(Cluster_specific_DEGs), file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_AT_ID.txt"))
      
      CS_DEGs_description = merge(Cluster_specific_DEGs, description, by.x = "gene_ID", by.y = "AT_ID", all.x = TRUE)[-c(8, 9)]
      CS_DEGs_description = CS_DEGs_description[order(CS_DEGs_description$pct.2, decreasing = FALSE), ]
      rownames(CS_DEGs_description) = CS_DEGs_description$gene_ID
      
      # Save the file
      write.csv(CS_DEGs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_description.csv"), row.names = FALSE)

      # Description - Cluster specific TFs 
      CS_TFs_description = CS_DEGs_description[!is.na(CS_DEGs_description$TF), ]
      
      # Save the file
      write.csv(CS_TFs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_description.csv"), row.names = FALSE)
      
      # Write the TF hirsuta ID file
      write_lines(CS_TFs_description$gene_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_AT_ID.txt"))
      
      AT_temp_cluster_DEG_list[[i]] <- rownames(Cluster_DEGs)
      names(AT_temp_cluster_DEG_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
      
      AT_specific_cluster_marker_list[[i]] <- rownames(Cluster_specific_DEGs)
      names(AT_specific_cluster_marker_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
      
    }
    
    else {
      # Prepare the file of cluster differentially expressed gene
      Cluster_DEGs = Cluster_DEGs[order(Cluster_DEGs[, colnames(temp_all_pct2)[1]], Cluster_DEGs[, colnames(temp_all_pct2)[2]]), ]
      Cluster_DEGs$gene_ID = rownames(Cluster_DEGs)
      
      # Save the file
      write.csv(Cluster_DEGs, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs.csv"), row.names = FALSE)
      
      # Write a text file containing Cardamine gene IDs
      write_lines(rownames(Cluster_DEGs), file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs_AT_ID.txt"))
      
      # Create ortho genes list and gene description file
      description = gene_description_file[rownames(gene_description_file) %in% rownames(Cluster_DEGs), -c(3, 4)]
      
      
      # Description - All DEGs identified for the cluster
      All_DEGs_description = merge(Cluster_DEGs, description, by.x = "gene_ID", by.y = "AT_ID", all.x = TRUE)
      All_DEGs_description = All_DEGs_description[order(All_DEGs_description[, colnames(temp_all_pct2)[1]], All_DEGs_description[, colnames(temp_all_pct2)[2]]), ]
      rownames(All_DEGs_description) = All_DEGs_description$gene_ID
      
      # Save the file
      write.csv(All_DEGs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_DEGs_description.csv"), row.names = FALSE)
      
      # Let's get all TFs information
      All_TFs_description = All_DEGs_description[!is.na(All_DEGs_description$TF), ]
      
      # Save the file
      write.csv(All_TFs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_TF_description.csv"), row.names = FALSE)
      
      # Write the TF hirsuta ID file
      write_lines(All_TFs_description$gene_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_all_TF_AT_ID.txt"))
      
      # Select cluster specific markers using the specific conserved marker finder function
      Cluster_specific_DEGs = specific_conserved_marker_finder(Cluster_DEGs, 
                                                               max.pct2.detection = 0.1, 
                                                               include.pct.difference = TRUE, 
                                                               pct.difference = 0.3, 
                                                               store.marker.file = FALSE)
      Cluster_specific_DEGs = Cluster_specific_DEGs[order(Cluster_specific_DEGs[, colnames(temp_all_pct2)[1]], Cluster_specific_DEGs[, colnames(temp_all_pct2)[2]]), ]
      Cluster_specific_DEGs$gene_ID = rownames(Cluster_specific_DEGs)
      
      # Save the file
      write.csv(Cluster_specific_DEGs, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers.csv"), row.names = FALSE)
      
      # Write a text file containing Cardamine gene IDs
      write_lines(rownames(Cluster_specific_DEGs), file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_AT_ID.txt"))
      
      CS_DEGs_description = merge(Cluster_specific_DEGs, description, by.x = "gene_ID", by.y = "AT_ID", all.x = TRUE)
      CS_DEGs_description = CS_DEGs_description[order(CS_DEGs_description[, colnames(temp_all_pct2)[1]], CS_DEGs_description[, colnames(temp_all_pct2)[2]]), ]
      rownames(CS_DEGs_description) = CS_DEGs_description$gene_ID
      
      # Save the file
      write.csv(CS_DEGs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_description.csv"), row.names = FALSE)

      # Description - Cluster specific TFs 
      CS_TFs_description = CS_DEGs_description[!is.na(CS_DEGs_description$TF), ]
      
      # Save the file
      write.csv(CS_TFs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_description.csv"), row.names = FALSE)
      
      # Write the TF hirsuta ID file
      write_lines(CS_TFs_description$gene_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_AT_ID.txt"))

      AT_temp_cluster_DEG_list[[i]] <- rownames(Cluster_DEGs)
      names(AT_temp_cluster_DEG_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
      
      AT_specific_cluster_marker_list[[i]] <- rownames(Cluster_specific_DEGs)
      names(AT_specific_cluster_marker_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
      
    }
  }
  
  if (return_markers_list){
    markerList = list(cluster_all_degs_AT = AT_temp_cluster_DEG_list,
                      cluster_specific_marker_AT = AT_specific_cluster_marker_list)
    
    return(markerList)
  }
  
}
