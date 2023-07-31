# Visualize the cells in reduced dimension plot and highlights the cells by cluster and splitting variable
cells_per_cluster_split_viz <- function(seuratObject, 
                                        store_dir = NULL, 
                                        store_folder = "Reduced_dimension_visualization",
                                        split_variable = "Species", 
                                        dimension_reduction_name = "umap", 
                                        dims_to_plot = c(1, 2)){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", toupper(dimension_reduction_name), "_plot_cells_per_cluster_split_by", "_", split_variable))){
    dir.create(str_c(store_dir, "/", store_folder, "/", toupper(dimension_reduction_name), "_plot_cells_per_cluster_split_by", "_", split_variable), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/", toupper(dimension_reduction_name), "_plot_cells_per_cluster_split_by", "_", split_variable, "/")
  
  # What are the cluster labels
  Clusters_labels = levels(Idents(seuratObject))
  
  # Go through each cluster and generate a dotplot
  for (i in (1:length(levels(Idents(seuratObject))))){
    p <- DimPlot(seuratObject, reduction = dimension_reduction_name, split.by = split_variable, pt.size = 0.5, label = T, label.size = 10, repel = T, cells.highlight = WhichCells(seuratObject, idents = Clusters_labels[i]), 
                 cols.highlight = c(RColorBrewer::brewer.pal(9, name = "Reds")[8], "#999999"), sizes.highlight = 1) + 
      xlab(str_c(toupper(dimension_reduction_name), dims_to_plot[1])) + 
      ylab(str_c(toupper(dimension_reduction_name), dims_to_plot[2])) + 
      theme(
        legend.position = "none",
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 24, face = "bold"),
        strip.text = element_text(size = 24, face = "bold.italic"))
    ggsave(filename = str_c(temp_dir, "Cells_in_cluster_", Clusters_labels[i], "_dims_", paste(dims_to_plot, collapse = "_"), ".png"), plot = p, width = 32, height = 18, dpi = 300)
  }
}
