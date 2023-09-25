GEP_generator <- function(input_data, 
                          file_name_pattern = NULL, 
                          reduction_key = "^inmf", 
                          store_data = TRUE,
                          store_GEPs = FALSE,
                          store_dir = NULL, 
                          store_folder = "Basis_genes"){
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  
  if (!is.character(input_data)){
    data.object <- input_data
  } else {
    if (file_test("-f", input_data)){
      data.object <- loadRData(input_data)
    }
    else {
      # Lets get the saved file names
      data.object = str_sort(list.files(path = input_data, pattern = file_name_pattern), numeric = TRUE)
    }
  }
  
  # Check if the input data is an object or path with multiple files
  
  if (!is.character(data.object)){
    
    # check if it is a seurat object
    
    if (tolower(class(data.object)) == "seurat"){
      
      temp = names(data.object@reductions)
      
      temp = temp[str_detect(temp, pattern = reduction_key)]
      
      
      temp_mat <- as.data.frame(data.object@reductions[[temp]]@feature.loadings)
      colnames(temp_mat) <- gsub(pattern = unique(gsub(pattern = "[^[:alpha:]]", replacement = "", x = colnames(temp_mat))), replacement = "GEP", x = colnames(temp_mat))
      temp_mat[str_c("K_", ncol(data.object@reductions[[temp]]@feature.loadings))] = parse_number(colnames(temp_mat)[apply(temp_mat, 1, which.max)])
      temp_mat[ ,str_c("K_", ncol(data.object@reductions[[temp]]@feature.loadings))] = as.factor(temp_mat[ ,str_c("K_", ncol(data.object@reductions[[temp]]@feature.loadings))])
      
      Basis = temp_mat[ ,str_c("K_", ncol(data.object@reductions[[temp]]@feature.loadings)), drop = FALSE]
      
      if (store_GEPs == TRUE){
        for (j in c(1:length(levels(Basis[ , 1])))){
          temp_genes = rownames(Basis[Basis[ , 1] == j, 1, drop = FALSE])
          fileGenerator(temp_genes, fileName = str_c(temp_dir, "GEP_", j, ".txt"))
        }
      }
      
    }
    
    else if (tolower(class(data.object)) == "liger") {
      
      temp_mat <- as.data.frame(t(data.object@W))
      colnames(temp_mat) <- str_c("GEP_", parse_number(colnames(temp_mat)))
      temp_mat[str_c("K_", ncol(t(data.object@W)))] = parse_number(colnames(temp_mat)[apply(temp_mat, 1, which.max)])
      temp_mat[ , str_c("K_", ncol(t(data.object@W)))] = as.factor(temp_mat[ , str_c("K_", ncol(t(data.object@W)))])
      
      Basis = temp_mat[ , str_c("K_", ncol(t(data.object@W))), drop = FALSE]
      
      if (store_GEPs == TRUE){
        for (j in c(1:length(levels(Basis[ , 1])))){
          temp_genes = rownames(Basis[Basis[ , 1] == j, 1, drop = FALSE])
          fileGenerator(temp_genes, fileName = str_c(temp_dir, "GEP_", j, ".txt"))
        }
      }
    }
  }
  
  else {
    
    if (length(data.object) == 1) {
      
      temp_file_name = list.files(path = input_data, pattern = file_name_pattern)
      
      temp_data_object = loadRData(str_c(input_data, "/", temp_file_name))
      
      if (tolower(class(temp_data_object)) == "seurat"){
        
        temp = names(temp_data_object@reductions)
        
        temp = temp[str_detect(temp, pattern = reduction_key)]
        
        temp_mat <- as.data.frame(temp_data_object@reductions[[temp]]@feature.loadings)
        colnames(temp_mat) <- gsub(pattern = unique(gsub(pattern = "[^[:alpha:]]", replacement = "", x = colnames(temp_mat))), replacement = "GEP", x = colnames(temp_mat))
        temp_mat[str_c("K_", ncol(temp_data_object@reductions[[temp]]@feature.loadings))] = parse_number(colnames(temp_mat)[apply(temp_mat, 1, which.max)])
        temp_mat[ ,str_c("K_", ncol(temp_data_object@reductions[[temp]]@feature.loadings))] = as.factor(temp_mat[ ,str_c("K_", ncol(temp_data_object@reductions[[temp]]@feature.loadings))])
        
        Basis = temp_mat[ ,str_c("K_", ncol(temp_data_object@reductions[[temp]]@feature.loadings)), drop = FALSE]
        
        if (store_GEPs == TRUE){
          for (j in c(1:length(levels(Basis[ , 1])))){
            temp_genes = rownames(Basis[Basis[ , 1] == j, 1, drop = FALSE])
            fileGenerator(temp_genes, fileName = str_c(temp_dir, "GEP_", j, ".txt"))
          }
        }
        
      }
      
      else if (tolower(class(temp_data_object)) == "liger") {
        
        temp_mat <- as.data.frame(t(temp_data_object@W))
        colnames(temp_mat) <- str_c("GEP_", parse_number(colnames(temp_mat)))
        temp_mat[str_c("K_", ncol(t(temp_data_object@W)))] = parse_number(colnames(temp_mat)[apply(temp_mat, 1, which.max)])
        temp_mat[ , str_c("K_", ncol(t(temp_data_object@W)))] = as.factor(temp_mat[ , str_c("K_", ncol(t(temp_data_object@W)))])
        
        Basis = temp_mat[ , str_c("K_", ncol(t(temp_data_object@W))), drop = FALSE]
        
        if (store_GEPs == TRUE){
          for (j in c(1:length(levels(Basis[ , 1])))){
            temp_genes = rownames(Basis[Basis[ , 1] == j, 1, drop = FALSE])
            fileGenerator(temp_genes, fileName = str_c(temp_dir, "GEP_", j, ".txt"))
          }
        }
      }
    }
    
    else {
      
      temp_file_names = list.files(path = input_data, pattern = file_name_pattern)
      
      # Creating an empty list to store the grouping results
      df_list <- list()
      
      for(i in c(1:length(temp_file_names))){
        
        file_name = str_c(input_data, "/", temp_file_names[i]) # File location and name
        
        temp_data_object = loadRData(file_name)
        
        if (tolower(class(temp_data_object)) == "seurat"){
          
          temp = names(temp_data_object@reductions)
          
          temp = temp[str_detect(temp, pattern = reduction_key)]
          
          temp_mat <- as.data.frame(temp_data_object@reductions[[temp]]@feature.loadings)
          colnames(temp_mat) <- gsub(pattern = unique(gsub(pattern = "[^[:alpha:]]", replacement = "", x = colnames(temp_mat))), replacement = "GEP", x = colnames(temp_mat))
          temp_mat[str_c("K_", ncol(temp_data_object@reductions[[temp]]@feature.loadings))] = parse_number(colnames(temp_mat)[apply(temp_mat, 1, which.max)])
          temp_mat[ ,str_c("K_", ncol(temp_data_object@reductions[[temp]]@feature.loadings))] = as.factor(temp_mat[ ,str_c("K_", ncol(temp_data_object@reductions[[temp]]@feature.loadings))])
          
          geneList = temp_mat[ , str_c("K_", ncol(temp_data_object@reductions[[temp]]@feature.loadings))]
          names(geneList) = rownames(temp_mat)
          df <- data.frame(geneList)
          colnames(df) <- str_c("K_", ncol(temp_data_object@reductions[[temp]]@feature.loadings))
          
          df_list[[i]] <- df
          
        }
        
        else if (tolower(class(temp_data_object)) == "liger") {
          
          temp_mat <- as.data.frame(t(temp_data_object@W))
          colnames(temp_mat) <- str_c("GEP_", parse_number(colnames(temp_mat)))
          temp_mat[str_c("K_", ncol(t(temp_data_object@W)))] = parse_number(colnames(temp_mat)[apply(temp_mat, 1, which.max)])
          temp_mat[ , str_c("K_", ncol(t(temp_data_object@W)))] = as.factor(temp_mat[ , str_c("K_", ncol(t(temp_data_object@W)))])
          
          geneList = temp_mat[ , str_c("K_", ncol(t(temp_data_object@W)))]
          names(geneList) = rownames(temp_mat)
          df <- data.frame(geneList)
          colnames(df) <- str_c("K_", ncol(t(temp_data_object@W)))
          
          df_list[[i]] <- df
          
        }
        
      }
      
      Basis <- do.call(cbind.data.frame, df_list)
    }
    
  }
  
  if (store_data == TRUE) {
    
    save(Basis, file = str_c(temp_dir, "Basis.RData"))
  }
  
  else {
    return(Basis)
  }
} 
