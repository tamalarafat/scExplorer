# scale the coefficient to have unit scale of the coefficients (proportion of GEP usage sum to 1) of the cells.
scale_coefficient <- function(seurat_object_dir, 
                              file_name_pattern, 
                              store_dir = NULL,
                              store_folder = "Seurat_objects_with_scaled_coefficient",
                              umap_n_neighbors = 30,
                              umap_seed_use = 2
                              ){
  
  # Lets get the saved file names
  seuratObjects = str_sort(list.files(path = seurat_object_dir, pattern = file_name_pattern), numeric = TRUE)
  
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
  
  for (i in c(1:length(seuratObjects))){
    # Creating the file name with path location
    file_name = str_c(seurat_object_dir, seuratObjects[i])
    
    # Loading the liger object and storing it to a temporary variable
    integrated.data = loadRData(file_name) # Load the RData file and assign it to a variable using the function loadRData
    
    # Dataset information
    integrated.data$Datasets <- integrated.data$orig.ident
    integrated.data$Datasets <- factor(integrated.data$Datasets, levels = c("C1", "C2", "C3", "C5", "O1", "O2", "O3", "O7"), labels = c("COL-1", "COL-2", "COL-3", "COL-5", "OX-1", "OX-2", "OX-3", "OX-7"))
    
    # Species information
    integrated.data$Species <- integrated.data$orig.ident
    integrated.data$Species <- substr(integrated.data$Species, 1, nchar(as.character(integrated.data$Species)) - 1)
    integrated.data$Species <- factor(integrated.data$Species, levels = c("C", "O"), labels = c("Thaliana", "Hirsuta"))
    
    integrated.data@reductions[["sinmf"]] <-integrated.data@reductions[["inmf"]]
    integrated.data@reductions[["sinmf"]]@cell.embeddings <- scale(integrated.data@reductions[["sinmf"]]@cell.embeddings, center = TRUE, scale = TRUE)
    
    integrated.data <- RunUMAP(integrated.data, reduction = "sinmf", dims = 1:dim(integrated.data@reductions[["sinmf"]])[2], n.components = 2, n.neighbors = umap_n_neighbors, seed.use = umap_seed_use)
    integrated.data <- RunTSNE(integrated.data, reduction = "sinmf", dims = 1:dim(integrated.data@reductions[["sinmf"]])[2], check_duplicates = FALSE, dim.embed = 2)
    
    integrated.data <- FindNeighbors(integrated.data, reduction = "sinmf", dims = 1:dim(integrated.data@reductions[["sinmf"]])[2])
    
    integrated.data <- FindClusters(integrated.data, resolution = seq(0.1, 1.2, 0.1), n.start = 50, n.iter = 50)
    
    save(integrated.data, file = str_c(temp_dir, "seurat_object_of_K", dim(integrated.data@reductions[["sinmf"]])[2], ".RData"))
    
  }
}
