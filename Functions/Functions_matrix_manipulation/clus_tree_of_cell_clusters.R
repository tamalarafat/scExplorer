# Tree of the cell clusters from Seurat clustering analysis of multiple resolution parameter
clus_tree_of_cell_clusters <- function(seuratObject, 
                          marker_file, 
                          query_pattern, 
                          gene_ID_column,
                          gene_name_column = NULL, 
                          store_dir = NULL, 
                          store_folder = "Tree_of_cell_clusters_graph_based_clustering"){
  # gene_ID_column <- specify the column, containing the gene ids, with column name or index to visualize the expression
  
  # Let's create a directory to store the GO annotation results for different factorization
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  # Let's create a directory to store the GO annotation results
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # Check if the folder exists (nested to previously created folder), if not create one to store the figures for the clustering using the coefficient matrix with K (any)
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "Cell_clusters_factorization_K_", ncol(seuratObject@reductions$inmf)))){
    dir.create(str_c(store_dir, "/", store_folder, "/", "Cell_clusters_factorization_K_", ncol(seuratObject@reductions$inmf)), showWarnings = TRUE, recursive = F, mode = "0777")
  }
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/", store_folder, "/", "Cell_clusters_factorization_K_", ncol(seuratObject@reductions$inmf), "/")
  
  p <- clustree(seuratObject, prefix = query_pattern, node_text_size = 8, node_size = 14) + 
    labs(colour = "Resolution\nvalue", edge_colour = "Count", edge_alpha = "Proportion") + 
    theme(legend.key.size = unit(2,"line"),
          legend.text = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 22, colour = "black", face = "bold"))
  ggsave(filename = str_c(temp_dir, "Tree_of_cell_clusters_K_", ncol(seuratObject@reductions$inmf), ".png"), plot = p, width = 48, height = 22, dpi = 300)
  
  # Check if the folder exists (nested to previously created folder), if not create one to store the figures for the clustering using the coefficient matrix with K (any)
  if (!dir.exists(str_c(temp_dir, "Markers_expression_across_cell_clusters"))){
    dir.create(str_c(temp_dir, "Markers_expression_across_cell_clusters"), showWarnings = TRUE, recursive = F, mode = "0777")
  }
  
  temp_exp_dir = str_c(temp_dir, "Markers_expression_across_cell_clusters", "/")
  
  if (is.null(gene_name_column)){
    # Set the rownames of the file with the gene IDs, so we can go over the gene ids
    rownames(marker_file) <- marker_file[ , gene_ID_column]
  } else {
    rownames(marker_file) <- marker_file[ , gene_name_column]
  }
  
  for (i in c(1:nrow(marker_file))){
    tryCatch({
      p <- clustree(seuratObject, prefix = query_pattern, node_text_size = 8, node_colour = marker_file[i, gene_ID_column], node_colour_aggr = "mean", node_size = 14) + 
        labs(colour = str_c(rownames(marker_file[i, ]) ,"\n", "avg. Expression"), edge_colour = "Count", edge_alpha = "Proportion") + 
        theme(legend.key.size = unit(2,"line"),
              legend.text = element_text(size = 20, face = "bold"),
              legend.title = element_text(size = 22, colour = "black", face = "bold"))
      
      ggsave(filename = str_c(temp_exp_dir, rownames(marker_file[i, ]), "_expression.png"), plot = p, width = 28, height = 22, dpi = 300)
    }, 
    error = function(e){print(str_c(rownames(marker_file[i, ]), " gene is not present in the integrated data."))}
    )
  }
}
