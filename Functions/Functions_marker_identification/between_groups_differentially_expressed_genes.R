# Identify differentially expressed genes between datasets in each cluster
between_groups_differentially_expressed_genes <- function(seuratObject, 
                                                         between_group_variable, 
                                                         store_dir = NULL, 
                                                         store_folder = "Differentially_expressed_genes", 
                                                         min_cell_for_DE = 10, 
                                                         DEGtest = "bimod", 
                                                         lfc = 0.25, 
                                                         Positive_exp = TRUE, 
                                                         min.pct = 0.1,
                                                         min.diff.pct = -Inf,
                                                         verbose = TRUE,
                                                         only.pos = FALSE,
                                                         max.cells.per.ident = Inf,
                                                         random.seed = 1,
                                                         latent.vars = NULL,
                                                         slot = "data",
                                                         counts = numeric(),
                                                         cells.1 = NULL,
                                                         cells.2 = NULL,
                                                         features = NULL
){
  
  # merge_group_markers tells the function to find markers of each grouping variable and merge them
  # For replicates analysis we want to find markers for each replicate and find conserved markers between the replicates.
  # But in species analysis we want to find conserved markers between the species. In case if a cluster is specific to a species we want to use FindMarkers function to find markers for the species
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "Differentially_expressed_genes_between_", between_group_variable, "_per_cluster"))){
    dir.create(str_c(store_dir, "/", store_folder, "/", "Differentially_expressed_genes_between_", between_group_variable, "_per_cluster"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "Differentially_expressed_genes_between_", between_group_variable, "_per_cluster", "/", "DEG_between_", between_group_variable, "_DEtest_", DEGtest))){
    dir.create(str_c(store_dir, "/", store_folder, "/", "Differentially_expressed_genes_between_", between_group_variable, "_per_cluster", "/", "DEG_between_", between_group_variable, "_DEtest_", DEGtest), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }

  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/", store_folder, "/", "Differentially_expressed_genes_between_", between_group_variable, "_per_cluster", "/", "DEG_between_", between_group_variable, "_DEtest_", DEGtest, "/")
  
  # Adding the initial identity / cell cluster information in 
  seuratObject[["temp_identity"]] <- Idents(seuratObject)
  
  # Let's combine the cluster identity with the tissue or experiment or species information
  seuratObject[["combined_identity"]] <- factor(str_c(Idents(seuratObject), seuratObject@meta.data[[between_group_variable]], sep = "_"))
  
  # Let's store the initial identity labels
  initial_idents = str_sort(levels(seuratObject$temp_identity), numeric = TRUE)
  
  # Lets assign the cluster identity to the combined identity
  Idents(seuratObject) <- "combined_identity"
  
  # Let's have the identity classes sorted
  Idents(seuratObject) <- factor(Idents(seuratObject), levels = str_sort(levels(Idents(seuratObject)), numeric = TRUE))
  
  # get the cluster labels of the newly created identity labels for the cells by combining the true (original cell clusters labels) identity and replicate or species or condition information
  combined_idents <- str_sort(levels(Idents(seuratObject)), numeric = TRUE)
  
  # How many cells are in each identity
  temp_cell_per_cluster <- data.frame(table(seuratObject@meta.data[["combined_identity"]]))
  
  # Using dplyr to rename the column names, sorting the cluster idents, finally using tibble to assign rownames with the column cluster
  temp_cell_per_cluster <- temp_cell_per_cluster %>% rename_("cluster" = "Var1", "count" = "Freq") 
  rownames(temp_cell_per_cluster) <- temp_cell_per_cluster$cluster
  temp_cell_per_cluster = temp_cell_per_cluster["count"]
  # Using this table we can check if any cluster has a minimum number of cells to be considered for performing differential expression analysis
  
  # Get the original cluster identity information that we can use to search for the groups in each cluster
  temp_pattern = str_sort(unique(gsub("(_).*", combined_idents, replacement = "\\1")), numeric = TRUE)
  
  if (length(initial_idents) != length(temp_pattern)) {print("Something went wrong while creating the combined labels")}
  
  # go through each element of the temp_pattern (holds thhe initial cluster idents with a seprator "_", 1_, 2_ etc) to find replicates/conditions in cluster ident
  for (i in c(1:length(temp_pattern))){
    
    # Create directory for each cluster
    dir.create(str_c(temp_dir, "/Cluster_", temp_pattern[i], "markers"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    
    groups_per_cluster = combined_idents[grep(pattern = str_c("^", temp_pattern[i]), x = combined_idents)]
    
    # If there is only one replicate in a cluster, we will not perform a differential analysis between the groups/replicates in that cluster
    if (length(groups_per_cluster) == 1){
      temp_str = c(gsub(pattern = "[^[:alnum:] ]", replacement = "", x = str_split(groups_per_cluster[1], pattern = "_")[[1]][2]), str_split(groups_per_cluster[1], pattern = "_")[[1]][1])
      temp_text = str_c("Data", " ", temp_str[1], " ", "is the only group or replicate in cluster", " ", temp_str[2], " ", "and all the cells in this cluster belong to this group. For", " ", temp_str[1], ", DE analysis between groups or replicates in cluster", " ", temp_str[2], " ", "can not be performed.")
      cat(temp_text, file = str_c(temp_dir, "/Cluster_", temp_pattern[i], "markers", "/No_markers_cell_C", temp_str[2], "_", temp_str[1], ".txt"))
    }
    
    else {
      for (j in c(1:length(groups_per_cluster))){
        # In the current cluster, if the number of cells of the current replicate is less than the defined thresold "min_cell_for_DE", we will not perfrom DE analysis for this replicate
        if (temp_cell_per_cluster[groups_per_cluster[j], ] < min_cell_for_DE){
          temp_str = c(gsub(pattern = "[^[:alnum:] ]", replacement = "", x = str_split(groups_per_cluster[j], pattern = "_")[[1]][2]), str_split(groups_per_cluster[j], pattern = "_")[[1]][1])
          temp_text = str_c("Data", " ", temp_str[1], " ", "in cluster", " ", temp_str[2], " ", "have less than", " ", min_cell_for_DE, " ", "cells. For", " ", temp_str[1], ", DE analysis between groups or replicates in cluster", " ", temp_str[2], " ", "was not performed.")
          cat(temp_text, file = str_c(temp_dir, "/Cluster_", temp_pattern[i], "markers", "/No_markers_cell_C", temp_str[2], "_", temp_str[1], ".txt"))
        }
        else {
          group_marker = FindMarkers(seuratObject, ident.1 = groups_per_cluster[j], ident.2 = groups_per_cluster[-j], test.use = DEGtest, min.pct = min.pct, logfc.threshold = lfc, only.pos = Positive_exp)
          temp_str = c(gsub(pattern = "[^[:alnum:] ]", replacement = "", x = str_split(groups_per_cluster[j], pattern = "_")[[1]][2]), str_split(groups_per_cluster[j], pattern = "_")[[1]][1])
          do.call("<-", list(str_c("Cluster_", temp_str[2], "_", temp_str[1]), group_marker))
          save(list = str_c("Cluster_", temp_str[2], "_", temp_str[1]), file = str_c(temp_dir, "/Cluster_", temp_pattern[i], "markers", "/Cluster_", temp_str[2], "_", temp_str[1], "_DEGs.RData"))
        }
      }
    }
  }
}
