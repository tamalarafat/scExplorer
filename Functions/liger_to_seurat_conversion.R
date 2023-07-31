# Convert all the liger objects to seurat objects
liger_to_seurat_conversion <- function(liger_object_dir, 
                                       file_name_pattern,
                                       scale_coefficient = TRUE, 
                                       store_dir = NULL, 
                                       store_folder = NULL){
  # liger_object_dir - Path of the saved Liger object files
  # file_name_pattern <- File name preffix to search in the directory
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # Lets get the saved file names
  temp_file_names = str_sort(list.files(path = liger_object_dir, pattern = file_name_pattern), numeric = TRUE)
  
  # Creating a integrated seurat object to store all the data from different liger objects
  
  # Loading the liger object and storing it to a temporary variable
  temp_liger_object = loadRData(str_c(liger_object_dir, "/", temp_file_names[length(temp_file_names)]))
  
  # Converting the liger object to seurat object
  temp_seurat_object = ligerToSeurat(temp_liger_object, nms = NULL)
  
  temp_seurat_object@reductions[["inmf"]] <- NULL
  temp_seurat_object@reductions[["tsne"]] <- NULL
  
  for(i in c(1:length(temp_file_names))){
    file_name = str_c(liger_object_dir, "/", temp_file_names[i]) # File location and name
    
    Liger_object = loadRData(file_name) # Load the RData file and assign it to a variable using the function loaadRData
    
    integrated.data <- ligerToSeurat(Liger_object, nms = NULL)
    
    integrated.data@reductions[["lumap"]] <- integrated.data@reductions[["tsne"]]
    integrated.data@reductions[["tsne"]] <- NULL
    
    integrated.data@reductions[["sinmf"]] <-integrated.data@reductions[["inmf"]]
    integrated.data@reductions[["sinmf"]]@cell.embeddings <- scale(integrated.data@reductions[["sinmf"]]@cell.embeddings, center = TRUE, scale = TRUE)
    
    if (scale_coefficient == TRUE) {
      
      integrated.data <- RunUMAP(integrated.data, reduction = "sinmf", dims = 1:dim(integrated.data@reductions[["sinmf"]])[2], n.components = 2)
      integrated.data <- RunTSNE(integrated.data, reduction = "sinmf", dims = 1:dim(integrated.data@reductions[["sinmf"]])[2], check_duplicates = FALSE, dim.embed = 2)
      
      integrated.data <- FindNeighbors(integrated.data, reduction = "sinmf", dims = 1:dim(integrated.data@reductions[["sinmf"]])[2])
      
      integrated.data <- FindClusters(integrated.data, resolution = seq(0.1, 1.2, 0.1), n.start = 50, n.iter = 50)
      
      save(integrated.data, file = str_c(temp_dir, "seurat_object_of_K_", dim(integrated.data@reductions[["sinmf"]])[2], ".RData"))
    }
    
    else {
      
      integrated.data <- RunUMAP(integrated.data, reduction = "inmf", dims = 1:dim(integrated.data@reductions[["inmf"]])[2], n.components = 2)
      integrated.data <- RunTSNE(integrated.data, reduction = "inmf", dims = 1:dim(integrated.data@reductions[["inmf"]])[2], check_duplicates = FALSE, dim.embed = 2)
      
      integrated.data <- FindNeighbors(integrated.data, reduction = "inmf", dims = 1:dim(integrated.data@reductions[["inmf"]])[2])
      
      integrated.data <- FindClusters(integrated.data, resolution = seq(0.1, 1.2, 0.1), n.start = 50, n.iter = 50)
      
      save(integrated.data, file = str_c(temp_dir, "seurat_object_of_K_", dim(integrated.data@reductions[["inmf"]])[2], ".RData"))
    }
    
    temp_seurat_object[[str_c("QNK", dim(integrated.data@reductions[["inmf"]])[2])]] <- Liger_object@clusters
    temp_seurat_object@reductions[[str_c("inmf", dim(integrated.data@reductions[["inmf"]])[2])]] <- integrated.data@reductions$inmf
    temp_seurat_object@reductions[[str_c("lumap", dim(integrated.data@reductions[["inmf"]])[2])]] <- integrated.data@reductions$LUMAP
  }
  
  # RUN -> the UMAP function
  temp_seurat_object <- RunUMAP(temp_seurat_object, reduction = str_c("inmf", parse_number(temp_file_names[length(temp_file_names)])), dims = 1:dim(temp_seurat_object@reductions[[str_c("inmf", parse_number(temp_file_names[length(temp_file_names)]))]])[2], n.components = 2)
  temp_seurat_object <- RunTSNE(temp_seurat_object, reduction = str_c("inmf", parse_number(temp_file_names[length(temp_file_names)])), dims = 1:dim(temp_seurat_object@reductions[[str_c("inmf", parse_number(temp_file_names[length(temp_file_names)]))]])[2], check_duplicates = FALSE, dim.embed = 2)
  
  temp_seurat_object <- FindNeighbors(temp_seurat_object, reduction = str_c("inmf", parse_number(temp_file_names[length(temp_file_names)])), dims = 1:dim(temp_seurat_object@reductions[[str_c("inmf", parse_number(temp_file_names[length(temp_file_names)]))]])[2])
  
  temp_seurat_object <- FindClusters(temp_seurat_object, resolution = seq(0.1, 1.2, 0.1), n.start = 50, n.iter = 50)
  
  integrated.data <- temp_seurat_object
  save(integrated.data, file = str_c(temp_dir, "seurat_object_with_all_factorized_K.RData"))
}