# For a series of seurat objects
calc_alignment_and_agreement <- function(objects_dir, 
                                         file_name_pattern = "Liger_object_K", 
                                         store_dir = NULL,
                                         store_folder = "Alignment_and_agreement_scores"
                                         ){
  
  # Lets get the saved file names
  list_objects = str_sort(list.files(path = objects_dir, pattern = file_name_pattern), numeric = TRUE)
  
  # Let's create a directory to store the GO annotation results for different factorization
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  # Let's create a directory to store the GO annotation results
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # Creating an empty list to store the grouping results
  df_list <- list()
  
  for (i in c(1:length(list_objects))){
    # Creating the file name with path location
    file_name = str_c(objects_dir, "/", list_objects[i])
    
    # Loading the liger object and storing it to a temporary variable
    Liger_object = loadRData(file_name) # Load the RData file and assign it to a variable using the function loadRData
    
    alignment_score = calcAlignment(Liger_object)
    
    agreement_score = calcAgreement(Liger_object, ndims = ncol(Liger_object@H.norm))
    
    df_list[[i]] = data.frame("Alignment" = alignment_score, "Agreement" = agreement_score)
    
    names(df_list)[i] <- str_c("GEP_", ncol(Liger_object@H.norm))
    
  }
  
  alignment_and_agreement <- do.call(rbind.data.frame, df_list)
  
  save(alignment_and_agreement, file = str_c(temp_dir, "alignment_and_agreement_score.RData"))
}
