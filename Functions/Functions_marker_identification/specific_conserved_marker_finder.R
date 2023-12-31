specific_conserved_marker_finder <- function(DEG_file,  
                                             max_pct2_detection = 0.1,
                                             include_pct_diff = FALSE,
                                             pct_diff = 0.5,
                                             file_name_pattern = "Cluster_",
                                             store_outputs = FALSE,
                                             store_dir = NULL, 
                                             store_folder = "Cluster_specific_conserved_marker",
                                             marker_file_name_preffix = "Specific_CM_"){
  
  # DEG_file = Differentially expressed genes identified previously using the "FindMarkers" function from seurat (identifies DEG for whole cluster); file can be a DEG dataframe or a list containing more than one DEG files
  # store_dir = path of the directory to store the files
  # store_folder = name of the folder to store the files
  # max_pct2_detection = primary criteria to select cluster specific marker genes
  # include_pct_diff = whether to consider pct_diff as an additional criteria to select cluster-specific mmarker genes
  # pct_diff = threshold of the difference in expression detection of a gene
  # max.pct.rest = additional criteria that can be jointly used to identify markers that are detected more than the defined max_pct2_detection threshold but highly expressed in the target cluster
  # store_outputs = whether to save the files or not
  
  # Lets convert the "DEG_file" to a list so that we can iterate over the list item in case of more than one DEG files
  
  # Creating necessary storing space to store the results
  
  if (store_outputs == TRUE) {
    
    if (missing(store_dir)){
      store_dir = getwd()
    }
    
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    temp_dir = str_c(store_dir, "/", store_folder, "/")
    
  }
  
  # Initializing an empty list to store all the DEG files
  temp_deg_list = list()
  
  # Check the files and load them
  if (class(DEG_file) == "character"){
    DEG_files = str_sort(list.files(DEG_file, pattern = file_name_pattern), numeric = TRUE)
    
    for (i in c(1:length(DEG_files))){
      
      CM = loadRData(str_c(DEG_dir, "/", DEG_files[i]))
      
      temp_deg_list[[i]] <- CM
      
      names(temp_deg_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]), "_CM")
    }
  }
  
  else if (class(DEG_file) != "character"){
    
    if (class(DEG_file) != "list"){
      temp_deg_list = list(DEG_file)
      
      if (length(temp_deg_list) == 1){
        
      names(temp_deg_list) <- "Cluster_specific_CM"
      }
    }
  }
  
  # Let's create a list to store the specific marker files to the list
  temp_list = list()
  
  # Get the names of the list items or DEG files
  temp_names = names(temp_deg_list)
  
  for (i in c(1:length(temp_names))){
    temp_deg_file = temp_deg_list[[temp_names[i]]]
    
    # max_pct2_detection = 0.1 (default); as we are looking to find specific conserved markers, there are more than one pct.2 columns (for each grouping variable (replicates / species / datasets)
    # The max_pct2_detection threshold identifies genes that are detected in a similar fashion across the grouping variable.
    # For example, if a gene is detected less than 10% cells in the other cells (pct.2) for each grouping varibale, the gene will be selected as specific conserved marker.

    # Extract all the columns in the DEG file containing pct.2 information
    temp_all_pct2 = temp_deg_file[ , str_detect(colnames(temp_deg_file), pattern = "pct.2"), drop = FALSE]
    
    if (ncol(temp_all_pct2) < 2){
      temp_specific_marker = specific_marker_finder(DEG_file_dir = temp_deg_file, store_outputs = FALSE)
      temp_specific_marker$gene_ID = rownames(temp_specific_marker)
      
      # checks if the user want to save the markers file in the directory or not
      if (store_outputs == TRUE){
        # Save the file
        write.csv(temp_specific_marker, file = str_c(temp_dir, "Cluster_", parse_number(temp_names[i]), "_specific_markers.csv"), row.names = TRUE)
      }
      
      else {
        return(temp_specific_marker)
      }
    }
    
    else {
      # get all the genes that are detected in user specified threshold, proportion of rest of the cells (vs target cluster) less than or equal to max_pct2_detection threshold.
      temp_primary_markers = temp_all_pct2[rowSums(temp_all_pct2 <= max_pct2_detection) == ncol(temp_all_pct2), ]
      
      # Get the row names == gene names; the primary list of specific markers
      temp_primary_gene_set = as.character(rownames(temp_primary_markers))
      
      # The chunk below is created to use the output (if we want to include pct difference criteria to look for markers) separately so that we can add the genes with the intended cluster-specific conserved markers list
      # include_pct_diff = FALSE; whether to include the pct difference criteria to look for cluster-specific conserved marker genes; if yes, use the two arguments below
      # pct_diff = 0.5;to include markers that has detection difference of atleast the specified threshold
      
      # Get all the columns in the DEG file that has pct information
      temp_pct_names = colnames(temp_deg_file)[str_detect(colnames(temp_deg_file), pattern = "pct")]
      
      # For each group we have pct.1, and pct.2 information. And we want to get the difference between these two proportion for each group.
      # Let's get the group names in the DEG file
      temp_group_names = unique(str_extract(temp_pct_names, pattern = "[^_]*"))
      
      # Let's create a list to store all the differnce outputs
      temp_diff_list = list()
      
      # Let's iterate over each group to calculate the difference and add them to the list as dataframe
      for (j in c(1:length(temp_group_names))){
        # Extract only the pct columns
        temp_object = temp_deg_file[ , c(str_c(temp_group_names[j], "_pct.1"), str_c(temp_group_names[j], "_pct.2"))]
        
        # Calculate the difference between the two pcts
        temp_object[[str_c("pct_diff_", temp_group_names[j])]] = temp_object[, str_c(temp_group_names[j], "_pct.1")] - temp_object[, str_c(temp_group_names[j], "_pct.2")]
        
        # Extract the pct difference column
        temp_diff = temp_object[, str_c("pct_diff_", temp_group_names[j]), drop = FALSE]
        
        # add the column to the empty list which will be used later create a dataframe with only the pct difference information
        temp_diff_list[[j]] <- temp_diff
      }
      
      # Let's transform the list object to a dataframe
      pct_diff_df = do.call(cbind.data.frame, temp_diff_list)
      
      # Keep only those genes that meet the criteria specified by the used with pct_diff threshold
      pct_diff_df <- pct_diff_df[rowSums(pct_diff_df >= pct_diff) == ncol(pct_diff_df), ]
      
      # Get the row names == gene names; get the additional gene list
      temp_additional_gene_set = as.character(rownames(pct_diff_df))
      
      # checks if we want to include the pct_diff criteria; include_pct_diff == TRUE will add the additional markers with the primary set of markers
      if (include_pct_diff == TRUE){
        temp_final_markers = unique(c(temp_primary_gene_set, temp_additional_gene_set))
        temp_specific_marker = temp_deg_file[temp_final_markers, ]
        temp_specific_marker$gene_ID = rownames(temp_specific_marker)
      }
      else {
        temp_specific_marker = temp_deg_file[temp_primary_gene_set, ]
        temp_specific_marker$gene_ID = rownames(temp_specific_marker)
      }
      
      # checks if the user want to save the markers file in the directory or not
      if (store_outputs == TRUE){
        write.csv(temp_specific_marker, file = str_c(temp_dir, "Cluster_", parse_number(temp_names[i]), "_specific_conserved_markers.csv"), row.names = TRUE)
      }
      
      else {
        return (temp_specific_marker)
      }
    }
  }
}

