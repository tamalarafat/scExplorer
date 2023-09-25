# 2 - Function to cluster the genes from the basis matrix for each factorization
cluster_basis_genes <- function(seurat_object_dir, 
                                reduction_name_pattern,
                                store_data = FALSE,
                                store_dir = NULL, 
                                store_folder = NULL){
  
  # liger_object_dir - Path of the saved Liger object files
  # file_name_pattern <- File name preffix to search in the directory
  
  # load the seurat object with all the factorized basis components included
  seurat_object = loadRData(seurat_object_dir)
  
  ### This function will generate a dataframe with cluster/coexpression list information of the genes
  temp = names(seurat_object@reductions)
  temp = str_sort(temp[str_detect(temp, pattern = reduction_name_pattern)], numeric = TRUE)
  
  # Creating an empty list to store the grouping results
  df_list <- list()
  
  for (i in c(1:length(temp))){
    temp_mat <- as.data.frame(seurat_object@reductions[[temp[i]]]@feature.loadings)
    temp_mat[str_c("LF_", parse_number(temp[i]))] = colnames(temp_mat)[apply(temp_mat, 1, which.max)]
    temp_mat[str_c("K_", parse_number(temp[i]))] = temp_mat[str_c("LF_", parse_number(temp[i]))]
    temp_mat[, str_c("K_", parse_number(temp[i]))] = parse_number(temp_mat[ ,str_c("K_", parse_number(temp[i]))])
    temp_mat["genes"] = rownames(temp_mat)
    temp_mat[ ,str_c("K_", parse_number(temp[i]))] = as.factor(temp_mat[ ,str_c("K_", parse_number(temp[i]))])
    
    geneList = temp_mat[ ,str_c("K_", parse_number(temp[i]))]
    names(geneList) = rownames(temp_mat)
    df <- data.frame(geneList)
    colnames(df) <- str_c("K_", parse_number(temp[i]))
    
    df_list[[i]] <- df
  }
  
  Basis <- do.call(cbind.data.frame, df_list)
  
  if (store_data == TRUE) {
    if (missing(store_dir)){
      store_dir = getwd()
    }
    
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # storing the directory information in a temporary variable
    temp_dir = str_c(store_dir, "/", store_folder, "/")
    
    save(Basis, file = str_c(temp_dir, "Basis.RData"))
  }
  
  else {
    return(Basis)
  }
}