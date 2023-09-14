# Get candidate markers from the identified differentially expressed genes for the clusters
candidate_markers_DEGs <- function(DEG_file,
                                   file_name_pattern = "Cluster",
                                   max_pct2_detection = 0.1,
                                   pct_diff = 0.3,
                                   include_pct_diff = TRUE,
                                   find_candidates = 10,
                                   cluster_ID = NULL,
                                   specify_gene_ID = NULL,
                                   incorporate_column_name = NULL,
                                   from_specified_list = FALSE,
                                   from_non_specified_list = FALSE,
                                   combine_categories = FALSE,
                                   store_outputs = TRUE,
                                   store_dir = NULL, 
                                   store_folder = "Candidate_markers_of_cluster_DEGs"){
  
  # Creating necessary storing space to store the results
  
  if (store_outputs == TRUE) {
    # Output directory
    if (missing(store_dir)){
      store_dir = getwd()
    }
    
    # Output folder
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Assigning the output directory path to a variable
    temp_dir = str_c(store_dir, "/", store_folder, "/")
  }
  
  
  # This part will take the argument input and get the genes of the GEP or GEPs and assign them into a list that will be used later to iterate over the GEP IDs.
  
  # A nested function to select candidates from the markers file
  select_candidates <- function(marker_file, number_of_candidates = find_candidates) {
    
    n_pct = colnames(marker_file)[str_detect(colnames(marker_file), pattern = "pct.2")]
    
    if (!"gene_ID" %in% colnames(marker_file)){
      marker_file$gene_ID = rownames(marker_file)
    }
    
    
    if (length(n_pct) == 1){
      # Let's define the candidates selection criteria
      marker_file = marker_file[order(marker_file[["pct.2"]], decreasing = FALSE), ]
      
      if (nrow(marker_file) == 0) {
        print("No TFs")
      } 
      
      else if (nrow(marker_file) <= number_of_candidates & nrow(marker_file) != 0){
        marker_file = marker_file
      } 
      
      else {
        if (nrow(marker_file[marker_file$avg_log2FC >= 1,]) < number_of_candidates & nrow(marker_file[marker_file$avg_log2FC >= 1,]) > 0){
          subset1 = marker_file[marker_file$avg_log2FC >= 1,]
          subset2 = marker_file[marker_file$avg_log2FC < 1,]
          subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
          
          if (nrow(subset2) <= (number_of_candidates - nrow(subset1))) {
            marker_file = rbind.data.frame(subset1, subset2)
          } else {
            subset2 = head(subset2, (number_of_candidates - nrow(subset1)))
            marker_file = rbind.data.frame(subset1, subset2)
          }
        }
        else if (nrow(marker_file[marker_file$avg_log2FC >= 1,]) == 0){
          marker_file = marker_file[order(marker_file$avg_log2FC, decreasing = TRUE), ]
          marker_file = marker_file[c(1:number_of_candidates), ]
        }
        else {
          marker_file = marker_file[marker_file$avg_log2FC >= 1, ]
          marker_file = marker_file[c(1:number_of_candidates), ]
        }
      }
      
    }
    
    else if (length(n_pct) == 2){
      marker_file = marker_file[order(marker_file[, n_pct[1]], marker_file[, n_pct[2]]), ]
      rownames(marker_file) = marker_file$gene_ID
      
      avg_FC_names = colnames(marker_file)[str_detect(colnames(marker_file), pattern = "avg_log2FC")]
      
      if (nrow(marker_file) == 0) {
        print("No TFs")} 
      
      else if (nrow(marker_file) <= number_of_candidates & nrow(marker_file) != 0){
        marker_file = marker_file
      } 
      
      else {
        if (nrow(marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]) < number_of_candidates & 
            nrow(marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]) > 0){
          subset1 = marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]
          subset2 = marker_file[marker_file[[avg_FC_names[1]]] < 1 & marker_file[[avg_FC_names[2]]] < 1, ]
          subset2 = subset2[order(subset2[, avg_FC_names[1]], subset2[, avg_FC_names[2]], decreasing = TRUE), ]
          
          if (nrow(subset2) <= (number_of_candidates - nrow(subset1))) {
            marker_file = rbind.data.frame(subset1, subset2)
          } 
          
          else {
            subset2 = head(subset2, (number_of_candidates - nrow(subset1)))
            marker_file = rbind.data.frame(subset1, subset2)
          }
        }
        else if (nrow(marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]) == 0){
          marker_file = marker_file[order(marker_file[, avg_FC_names[1]], marker_file[, avg_FC_names[2]], decreasing = TRUE), ]
          marker_file = marker_file[c(1:number_of_candidates), ]
        }
        
        else {
          marker_file = marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]
          marker_file = marker_file[c(1:number_of_candidates), ]
        }
      }
    }
    
    # Add the source of the markers
    marker_file$source = "DEGs"
    
    # Return the marker file
    return(marker_file)
  }
  
  # Create an empty list to store the deg files
  temp_deg_list = list()
  
  if (is.character(DEG_file)) {
    
    if (!file_test("-f", DEG_file)) {
      
      # get the file names
      file_names = str_sort(list.files(DEG_file, pattern = file_name_pattern), numeric = TRUE)
      
      for (i in c(1:length(file_names))){
        
        marker_file = loadRData(str_c(DEG_file, "/", file_names[i]))
        
        temp_deg_list[[i]] <- marker_file
        
        if (grepl("[[:digit:]]", file_names[i]) != TRUE){
          names(temp_deg_list)[i] <- sub("(.*)\\.RData$", "\\1", file_names[i])
        }
        
        else {
          names(temp_deg_list)[i] <- str_c("cluster_", parse_number(file_names[i]), "_degs")
        }
      }
    }
    
    else {
      
      # If the input is a caharacter vector containing the file name (the file has to be ".RData" format)
      marker_file = loadRData(DEG_file)
      
      # Add the file to the list (containing the degs table)
      temp_deg_list[[1]] <- marker_file
      
      # Get the name of the deg file
      file_names = basename(DEG_file)
      
      # Check if the file name contains any digit, indicator of a cluster number. If yes, get the digit, else, use the file name
      
      # if no digit, get the file name
      if (grepl("[[:digit:]]", file_names) != TRUE){
        names(temp_deg_list)[1] <- sub("(.*)\\.RData$", "\\1", file_names)
      }
      
      # else use the digit to name the list item
      else {
        names(temp_deg_list)[1] <- str_c("cluster_", parse_number(file_names), "_degs")
      }
    }
  } 
  
  # If the deg file is not a character, check if it's a list or data frame
  else {
    
    # If the input is not a list but a data.frame
    if (!is.list(class(DEG_file))){
      temp_deg_list[[1]] <- DEG_file
      
      # Assign the name to the first item of the
      names(temp_deg_list)[1] <- ifelse(!is.null(cluster_ID), str_c("cluster_", cluster_ID, "_degs"), "cluster_degs")
    }
    
    # else, assign the input deg list to the temporary variable
    else {
      temp_deg_list = DEG_file
      
      # If the names of the list items is null set the names
      if (is.null(temp_deg_list)) {
        names(temp_deg_list) <- str_c("list_item_", seq(1:length(temp_deg_list)))
      }
    }
  }
  
  
  # create an empty list to store the candidate markers files
  temp_list = list()
  
  
  for (i in c(1:length(temp_deg_list))) {
    
    # Get the deg table and assign it to the variable
    cluster_marker = temp_deg_list[[i]]
    
    # Check if the deg table obtained from integrated/comparative analysis (conserved) or single species analysis (typical).
    
    # Get all the "pct.2" from the deg table 
    n_pct = colnames(cluster_marker)[str_detect(colnames(cluster_marker), pattern = "pct.2")]
    
    # Get cluster specific markers for the cell clusters
    if (length(n_pct) == 1){
      # Select specific marker from the deg file - typical deg table for a single dataset
      cluster_marker = specific_marker_finder(DEG_file_dir = cluster_marker, max_pct2_detection = max_pct2_detection, pct_diff = pct_diff, include_pct_diff = include_pct_diff, store_outputs = FALSE)
    }
    
    else if (length(n_pct) == 2) {
      # Select specific marker from the deg file - conserved
      cluster_marker = specific_conserved_marker_finder(DEG_file = cluster_marker, max_pct2_detection = max_pct2_detection, pct_diff = pct_diff, include_pct_diff = include_pct_diff, store_outputs = FALSE)
    }
    
    
    ### This part checks if the user specified a set of genes from which they want to select candidate markers
    
    # Check if user specified a set of genes to look for candidates from the DEGs file
    if (!is.null(specify_gene_ID)) {
      
      ## This part for the specified genes set
      
      # First, subset the marker file with the specified gene IDs
      specified_cluster_marker = cluster_marker[cluster_marker$gene_ID %in% specify_gene_ID, ]
      
      if (nrow(specified_cluster_marker) == 0) {
        specified_cluster_candidates = specified_cluster_marker
      } 
      
      else {
        # Find candidate markers for the specified genes
        specified_cluster_candidates = select_candidates(marker_file = specified_cluster_marker, number_of_candidates = find_candidates)
        
        # Check if user wants to add the information of the category to the table. If not missing add the column with the given name.
        if (!is.null(incorporate_column_name)) {
          specified_cluster_candidates[[incorporate_column_name]] <- "Yes"
        }
      }
      
      
      
      ## This part for the remaining genes
      
      # First, subset the marker file with the specified gene IDs
      non_specified_cluster_marker = cluster_marker[!cluster_marker$gene_ID %in% specify_gene_ID, ]
      
      if (nrow(non_specified_cluster_marker) == 0) {
        non_specified_cluster_candidates = non_specified_cluster_marker
      } 
      
      else {
        # Find candidate markers for the specified genes
        non_specified_cluster_candidates = select_candidates(marker_file = non_specified_cluster_marker, number_of_candidates = find_candidates)
        
        # Check if user wants to add the information of the category to the table. If not missing add the column with the given name.
        if (!is.null(incorporate_column_name)) {
          non_specified_cluster_candidates[[incorporate_column_name]] <- "No"
        }
      }
      
      
      ## Now check the requirements of the user.
      
      # If user wants to combine both categories
      if (combine_categories == TRUE) {
        
        # bind the rows of the two tables with row bind
        candidate_markers = rbind.data.frame(specified_cluster_candidates, non_specified_cluster_candidates)

      }
      
      # If user wants to find candidates only from the specified list of genes
      else if (from_specified_list == TRUE) {
        
        # Add the candidates data table to the list
        candidate_markers <- specified_cluster_candidates
        
      }
      
      # If user wants to find candidates only from the specified list of genes
      else if (from_non_specified_list == TRUE) {
        
        # Add the candidates data table to the list
        candidate_markers <- non_specified_cluster_marker

      }
    }
    
    else {
      
      # Find candidate markers for the specified genes
      candidate_markers = select_candidates(marker_file = cluster_marker, number_of_candidates = find_candidates)
      
    }
    
    
    # Rearrange the rows of the data table based on the pct.2 criteria. This will be used to rank the genes (lowest pct.2 gets the highest rank)
    
    # If the table is a typical deg table of a single species or dataset
    if (length(n_pct) == 1){
      candidate_markers = candidate_markers[order(candidate_markers[["pct.2"]], decreasing = FALSE), ]
    }
    
    # If the table is a conserved deg table of a two species or conditions
    else if (length(n_pct) == 2) {
      candidate_markers = candidate_markers[order(candidate_markers[, n_pct[1]], candidate_markers[, n_pct[2]]), ]
    }
    
    candidate_markers$rank = c(1:nrow(candidate_markers))
    
    # Add the candidates data table to the list
    temp_list[[i]] <- candidate_markers
    
    # Assign name to the list item
    names(temp_list)[i] <- names(temp_deg_list)[i]
  }
  
  # Check if the user provided any cluster ID
  if (!is.null(cluster_ID)){

    # Subset the list with the cluster ID
    temp_list = temp_list[grep(pattern = str_c("^", cluster_ID, "$", sep = ""), parse_number(names(temp_list)))]
  }

  # Check if store outputs is true. If TRUE store the files
  if (store_outputs == TRUE) {

    # Write csv files for each list item
    invisible(lapply(names(temp_list), function (i) {write.csv(temp_list[[i]], file = str_c(temp_dir, i, "_candidates.csv"), row.names = FALSE)}))
  }

  else {

    # If the list has only one item, return the data table
    ifelse(length(temp_list) == 1, return(temp_list[[1]]), return(temp_list))
  }
}

