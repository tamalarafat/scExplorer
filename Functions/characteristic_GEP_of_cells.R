# Get the significant GEP for a cluster, or a set of cells, or cells expressing a given gene
characteristic_GEP_of_cells <- function(seurat_object, 
                               target_ID, 
                               reduction_name = "^inmf", 
                               store_output = TRUE, 
                               cell_ids, 
                               gene_ID,
                               store_dir = NULL, 
                               store_folder = "Defining_GEP_of_cluster") {
  
  # Creating necessary storing space to store the results
  
  # Output directory
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  # Output folder
  
  if(!missing(target_ID)) {
    
    if (!dir.exists(str_c(store_dir, "/", store_folder, "_", target_ID))){
      dir.create(str_c(store_dir, "/", store_folder, "_", target_ID), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Assigning the output directory path to a variable
    temp_dir = str_c(store_dir, "/", store_folder, "_", target_ID, "/")
  }
  
  
  # Output folder
  
  if(missing(target_ID)) {
    
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Assigning the output directory path to a variable
    temp_dir = str_c(store_dir, "/", store_folder, "/")
  }
  
  
  # Check if user specified a cluster to check GEP usage
  if (!missing(target_ID)) {
    
  # Get the cell ids belonging to the cluster
  cluster_cells = WhichCells(seurat_object, idents = target_ID)
  
  # We dont need to load the liger object as the seurat object contains the normalized matrix
  inmf_mat = seurat_object@reductions[[reduction_name]]@cell.embeddings
  
  # Add colnames to the dataframe - Coefficient matrix
  colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))
  
  # Subset the coefficient matrix with the cluster cells
  inmf_mat = inmf_mat[cluster_cells, ]
  
  # Order the coefficient matrix to determine which GEP or GEPs the majority of the cluster cells are using.
  df_H <- order_factorized_matrix(inmf_mat)
  
  # For each cell, identify the GEP that is maximally utilized by the cell
  df_H$GEP <- apply(df_H, 1, function (x) parse_number(colnames(df_H)[which.max(x)]))
  
  # Generate a data frame that includes the total number of cells along with their corresponding maximum GEP.
  cluster_GEP_usage = data.frame(table(df_H$GEP))
  
  # Add colnames to the dataframe
  colnames(cluster_GEP_usage) = c("GEP", "cell_count")
  
  # Create an empty list to store the GEP usage by the cluster cells, provided set of cells, and cells expressing a given gene
  GEP_usage = list()
  
  # Store the GEP usage in the list
  GEP_usage[[1]] <- cluster_GEP_usage
  
  
  # Check max GEP for a provided set of cells
  if (!missing(cell_ids)) {
    
    # Get the cell ids present in the cluster
    provied_cells_in_cluster = intersect(cell_ids, cluster_cells)
    
    # subset the dataframe with only the provided cell ids
    cells_GEP = df_H[provied_cells_in_cluster, ]
    
    # Generate a data frame that includes the total number of cells along with their corresponding maximum GEP.
    cells_GEP = data.frame(table(cells_GEP$GEP))
    
    # Add colnames to the dataframe
    colnames(cells_GEP) = c("GEP", "cell_count")
    
    # Add the GEP usage to the list
    GEP_usage[[2]] <- cells_GEP
    
    # Check if the user wants to store the output or not
    if (store_output == TRUE) {
      
      # Store the output in the desired directory
      save(cells_GEP, file = str_c(temp_dir, "cluster_", target_ID, "_provied_cells_GEP_usage.RData"))
    }
    
    else {
      # return the GEP usage file
      return (GEP_usage)
    }
  }
  
  else if (!missing(gene_ID)) {
    
    # Check the STM expressing cells
    cells_with_expression_detection = names(GetAssayData(seurat_object, assay = "RNA", slot = "counts")[gene_ID, GetAssayData(seurat_object, assay = "RNA", slot = "counts")[gene_ID, ] != 0])
    
    # Get the cell ids present in the cluster
    cluster_expressing_cells = intersect(cells_with_expression_detection, cluster_cells)
    
    # subset the dataframe with only the provided cell ids
    cells_GEP = df_H[cluster_expressing_cells,]
    
    # Generate a data frame that includes the total number of cells along with their corresponding maximum GEP.
    cells_GEP = data.frame(table(cells_GEP$GEP))
    
    # Add colnames to the dataframe
    colnames(cells_GEP) = c("GEP", "cell_count")
    
    # Add the GEP usage to the list
    GEP_usage[[2]] <- cells_GEP
    
    # Check if the user wants to store the output or not
    if (store_output == TRUE) {
      
      # Store the output in the desired directory
      save(cells_GEP, file = str_c(temp_dir, "cluster_", target_ID, "_expressing_cells_GEP_usage.RData"))
    }
    
    else {
        
        # return the GEP usage file
        return (GEP_usage)
      }
  }
  
  else {
    
    # Check if the user wants to store the output or not
    if (store_output == TRUE) {
      
      # Store the output in the desired directory
      save(cluster_GEP_usage, file = str_c(temp_dir, "cluster_", target_ID, "_GEP_usage.RData"))
    }
    
    else {
      
      # return the output
      return (cluster_GEP_usage)
    }
  }
  }
  
  
  # If user provied cell ids or gene ID but not any cluster ID
  else {
    
    # Extract the coefficient matrix from the Seurat object
    inmf_mat = seurat_object@reductions[[reduction_name]]@cell.embeddings
    
    # Add colnames to the dataframe - Coefficient matrix
    colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))
    
    # Check max GEP for a provided set of cells
    if (!missing(cell_ids)) {
      
      # Subset the coefficient matrix with the cluster cells
      inmf_mat = inmf_mat[cell_ids, ]
      
      # Order the coefficient matrix to determine which GEP or GEPs the majority of the cluster cells are using.
      df_H <- order_factorized_matrix(inmf_mat)
      
      # For each cell, identify the GEP that is maximally utilized by the cell
      df_H$GEP <- apply(df_H, 1, function(x) parse_number(colnames(df_H)[which.max(x)]))
      
      # Generate a data frame that includes the total number of cells along with their corresponding maximum GEP.
      cells_GEP_usage = data.frame(table(df_H$GEP))
      
      # Add colnames to the dataframe
      colnames(cells_GEP_usage) = c("GEP", "cell_count")
      
      # Check if the user wants to store the output or not
      if (store_output == TRUE) {
        
        # Store the output in the desired directory
        save(cells_GEP_usage, file = str_c(temp_dir, "Provied_cells_GEP_usage.RData"))
      }
      
      else {
        
        # return the GEP usage file
        return (cells_GEP_usage)
      }
    }
    
    else if (!missing(gene_ID)) {
      
      # get the expressing cells
      cells_with_expression_detection = names(GetAssayData(seurat_object, assay = "RNA", slot = "counts")[gene_ID, GetAssayData(seurat_object, assay = "RNA", slot = "counts")[gene_ID, ] != 0])
      
      # Subset the coefficient matrix with the cluster cells
      inmf_mat = inmf_mat[cells_with_expression_detection, ]
      
      # Order the coefficient matrix to determine which GEP or GEPs the majority of the cluster cells are using.
      df_H <- order_factorized_matrix(inmf_mat)
      
      # For each cell, identify the GEP that is maximally utilized by the cell
      df_H$GEP <- apply(df_H, 1, function(x) parse_number(colnames(df_H)[which.max(x)]))
      
      # Generate a data frame that includes the total number of cells along with their corresponding maximum GEP.
      cells_GEP_usage = data.frame(table(df_H$GEP))
      
      # Add colnames to the dataframe
      colnames(cells_GEP_usage) = c("GEP", "cell_count")
      
      # Check if the user wants to store the output or not
      if (store_output == TRUE) {
        
        # Store the output in the desired directory
        save(cells_GEP_usage, file = str_c(temp_dir, gene_ID, "_expressing_cells_GEP_usage.RData"))
      }
      
      else {
        
        # return the GEP usage file
        return (cells_GEP_usage)
      }
    }
  }
}

