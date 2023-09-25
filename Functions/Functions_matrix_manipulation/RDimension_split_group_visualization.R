# Dimension reduction plot for the cells with species or replicate or dataset information
RDimension_split_group_visualization <- function(seuratObject,
                                                 store_dir = NULL,
                                                 store_folder = "Reduced_dimension_visualization",
                                                 dimension_reduction_name = "umap",
                                                 dims_to_plot = c(1, 2),
                                                 split_variable = "Species",
                                                 plot_axis_lines = TRUE){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", toupper(dimension_reduction_name), "_plot_split_by", "_", split_variable))){
    dir.create(str_c(store_dir, "/", store_folder, "/", toupper(dimension_reduction_name), "_plot_split_by", "_", split_variable), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/", toupper(dimension_reduction_name), "_plot_split_by", "_", split_variable, "/")
  
  Idents(seuratObject) <- split_variable
  
  if (plot_axis_lines == TRUE){
    p <- DimPlot(seuratObject, reduction = dimension_reduction_name, dims = dims_to_plot, label = T, split.by = split_variable, pt.size = 0.5, repel = T) + 
      NoLegend() + 
      theme(
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"), 
        axis.ticks.length = unit(.30, "cm"), 
        axis.text = element_text(size = 24, face = "bold"), 
        strip.text = element_text(size = 24, face = "bold.italic")) + 
      guides(colour = guide_legend(override.aes = list(size = 8)))
    
    ggsave(plot = p, filename = str_c(temp_dir, "Cells_", toupper(dimension_reduction_name), "_colored_with_", split_variable, "_dims_", paste(dims_to_plot, collapse = "_"), ".png"), width = 32, height = 18, dpi = 300)
  }
  
  else {
    p <- DimPlot(seuratObject, reduction = dimension_reduction_name, dims = dims_to_plot, label = T, split.by = split_variable, pt.size = 0.5, repel = T) + 
      NoLegend() + 
      theme(
        line = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        strip.text = element_text(size = 24, face = "bold.italic")) + 
      guides(colour = guide_legend(override.aes = list(size = 8)))
    
    ggsave(plot = p, filename = str_c(temp_dir, "Cells_", toupper(dimension_reduction_name), "_colored_with_", split_variable, "_dims_", paste(dims_to_plot, collapse = "_"), ".png"), width = 32, height = 18, dpi = 300)
  }
}
