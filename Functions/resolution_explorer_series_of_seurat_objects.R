# For a series of seurat objects
resolution_explorer_series_of_seurat_objects <- function(seurat_object_dir, 
                                                         file_name_pattern, 
                                                         store_dir = NULL,
                                                         store_folder = "Cells_genes_clusters",
                                                         ident_levels = NULL,
                                                         add_replicate_labels = NULL,
                                                         add_species_labels = NULL,
                                                         add_tissue_labels = NULL
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
  temp_main_dir = str_c(store_dir, "/", store_folder, "/")
  
  for (i in c(1:length(seuratObjects))){
    # Creating the file name with path location
    file_name = str_c(seurat_object_dir, "/", seuratObjects[i])
    
    # Loading the liger object and storing it to a temporary variable
    integrated.data = loadRData(file_name) # Load the RData file and assign it to a variable using the function loadRData
    
    # Dataset information
    integrated.data$Replicates <- integrated.data$orig.ident
    integrated.data$Replicates <- factor(integrated.data$Replicates, levels = ident_levels, labels = add_replicate_labels)
    
    # Species information
    integrated.data$Species <- integrated.data$orig.ident
    integrated.data$Species <- substr(integrated.data$Species, 1, nchar(as.character(integrated.data$Species)) - 1)
    integrated.data$Species <- factor(integrated.data$Species, levels = c("O", "S"), labels = c("Hirsuta", "AT-SAM"))

    # Tissue information
    integrated.data$Tissue <- integrated.data$orig.ident
    integrated.data$Tissue <- substr(integrated.data$Tissue, 1, nchar(as.character(integrated.data$Tissue)) - 1)
    integrated.data$Tissue <- factor(integrated.data$Tissue, levels = c("O", "S"), labels = c("Leaf", "Apex"))

    save(integrated.data, file = str_c(seurat_object_dir, "seurat_object_of_K_", dim(integrated.data@reductions[["inmf"]])[2], ".RData"))

    # Create a directory to store the results of the particular factorization
    # Check if the folder exists (nested to previously created folder), if not create one to store the figures for the clustering using the coefficient matrix with K (any)
    # if (!dir.exists(str_c(temp_main_dir, "Factor_K_", ncol(integrated.data@reductions$inmf)))){
    #   dir.create(str_c(temp_main_dir, "Factor_K_", ncol(integrated.data@reductions$inmf)), showWarnings = TRUE, recursive = F, mode = "0777")
    # }
    
    # temp_factor_dir = str_c(temp_main_dir, "Factor_K_", ncol(integrated.data@reductions$inmf), "/")
    
    resolution_explorer(seuratObject = integrated.data, store_dir = temp_main_dir, store_folder = str_c("Factor_K_", ncol(integrated.data@reductions$inmf)))
    
  }
}
