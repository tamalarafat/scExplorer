# Factorize to lower Ks
optimize_to_new_k <- function(optimized_liger_object, 
                              choice_of_new_K, 
                              store_dir = NULL, 
                              store_folder = NULL){
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # Here we load the LigerObject from the provided path and assign it to a temporary variable
  LigerObject = loadRData(optimized_liger_object)
  
  # Creating an empty list to store the quantile normalization results
  QNKs_Clusters = list()
  
  # In the for loop we will go through each K and factorize to new K
  for (i in c(1:length(choice_of_new_K))){
    Liger_object <- optimizeNewK(LigerObject, k.new = choice_of_new_K[i])
    
    Liger_object <- quantile_norm(Liger_object) # Perform quantile normalization
    
    Liger_object <- runUMAP(Liger_object)
    
    save(Liger_object, file = str_c(temp_dir, "Liger_object_K_", choice_of_new_K[i], ".RData"))
    
    temp = Liger_object@clusters
    
    QNKs_Clusters[[i]] = temp
    
    names(QNKs_Clusters)[i] = str_c("QNK", choice_of_new_K[i])
  }
  
  Liger_object = LigerObject
  
  save(Liger_object, file = str_c(temp_dir, "Liger_object_K_", dim(Liger_object@W)[1], ".RData"))
  
  temp = LigerObject@clusters
  
  qn_Clustering = data.frame(QNKs_Clusters)
  
  qn_Clustering = cbind.data.frame(qn_Clustering, temp)
  
  colnames(qn_Clustering)[dim(qn_Clustering)[2]] = str_c("QNK", dim(LigerObject@W)[1])
  
  save(qn_Clustering, file = str_c(temp_dir, "qn_Clustering.RData"))
}