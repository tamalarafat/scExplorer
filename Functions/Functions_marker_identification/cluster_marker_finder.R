# Identify differentially expressed genes for each cluster
cluster_marker_finder <- function(seuratObject, 
                                  store_dir = NULL, 
                                  store_folder = "Differentially_expressed_genes",
                                  min_cell_for_DE = 10, 
                                  DEGtest = "bimod", 
                                  lfc = 0.25, 
                                  positive_exp = TRUE, 
                                  min.pct = 0.1,
                                  min.diff.pct = -Inf,
                                  verbose = TRUE,
                                  max.cells.per.ident = Inf,
                                  random.seed = 1,
                                  latent.vars = NULL,
                                  slot = "data",
                                  counts = numeric(),
                                  cells.1 = NULL,
                                  cells.2 = NULL,
                                  features = NULL
){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "DEGs_of_each_cluster"))){
    dir.create(str_c(store_dir, "/", store_folder, "/", "DEGs_of_each_cluster"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "DEGs_of_each_cluster", "/", "DEG_DEtest_", DEGtest))){
    dir.create(str_c(store_dir, "/", store_folder, "/", "DEGs_of_each_cluster", "/", "DEG_DEtest_", DEGtest), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/", store_folder, "/", "DEGs_of_each_cluster", "/", "DEG_DEtest_", DEGtest, "/")
  
  # Let's store the initial identity labels
  initial_idents = str_sort(levels(Idents(seuratObject)), numeric = TRUE)
  
  # Let's have the identity classes sorted
  Idents(seuratObject) <- factor(Idents(seuratObject), levels = str_sort(levels(Idents(seuratObject)), numeric = TRUE))
  
  # How many cells are in each identity
  temp_cell_per_cluster <- data.frame(table(Idents(seuratObject)))
  
  # Using dplyr to rename the column names, sorting the cluster idents, finally using tibble to assign rownames with the column cluster
  temp_cell_per_cluster <- temp_cell_per_cluster %>% rename_("cluster" = "Var1", "count" = "Freq") 
  rownames(temp_cell_per_cluster) <- temp_cell_per_cluster$cluster
  temp_cell_per_cluster = temp_cell_per_cluster["count"]
  
  # Using this table we can check if any cluster has a minimum number of cells to be considered for performing differential expression analysis
  
  if (length(initial_idents) != nrow(temp_cell_per_cluster)) {print("Something went wrong with the cluster labels.")}
  
  # go through each element of the cluster ident
  for (i in c(1:length(initial_idents))){
    
    # In the current cluster, if the number of cells of the current replicate is less than the defined thresold "min_cell_for_DE", we will not perfrom DE analysis for this replicate
    if (temp_cell_per_cluster[initial_idents[i], ] < min_cell_for_DE){
      temp_text = str_c("Cluster", " ", initial_idents[i], " ", "have less than", " ", min_cell_for_DE, " ", "cells. For", " ", initial_idents[i], ", DE analysis was not performed.")
      cat(temp_text, file = str_c(temp_dir, "No_markers_cell_C", initial_idents[i], ".txt"))
    }
    else {
      cluster_marker = FindMarkers(seuratObject, ident.1 = initial_idents[i], test.use = DEGtest, min.pct = min.pct, logfc.threshold = lfc, only.pos = positive_exp)
      do.call("<-", list(str_c("Cluster_", initial_idents[i], "_DEGs"), cluster_marker))
      save(list = str_c("Cluster_", initial_idents[i], "_DEGs"), file = str_c(temp_dir, "Cluster_", initial_idents[i], "_DEGs.RData"))
    }
  }
}
