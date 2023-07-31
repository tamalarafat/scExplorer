# Function to cluster the genes from the basis matrix for each factorization
gene_clustering_of_nmfObject <- function(dir_path, fileName_pattern, store_output = FALSE){
  
  ### This function will generate a dataframe with cluster information of the genes
  nmfFiles = str_sort(list.files(path = nmfObject_dir, pattern = fileName_pattern), numeric = TRUE)
  
  # Creating an empty list to store the clustering result for each factorization
  df_list <- list()
  
  for (i in c(1:length(nmfFiles))){
    
    temp_file = str_c(dir_path, nmfFiles[i])
    
    classifier <- loadRData(temp_file)
    
    genes_cluster = predict(classifier, what = "rows")
    
    df <- data.frame(genes_cluster)
    
    colnames(df) <- str_c("K_", parse_number(nmfFiles[i]))
    
    df_list[[i]] <- df
    
    names(df_list)[i] <- str_c("K_", parse_number(nmfFiles[i]))
  }
  
  Basis <- do.call(cbind.data.frame, df_list)
  
  if (store_output) {
    save(Basis, file = str_c(dir_path, "Basis.RData"))
  }
  else {
    return(Basis)
  }
}