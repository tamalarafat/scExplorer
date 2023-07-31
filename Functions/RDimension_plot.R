# Dimension reduction plot for the cell clusters labels
RDimension_plot <- function(seuratObject, 
                            store_dir = NULL, 
                            store_folder = "Reduced_dimension_visualization",
                            dimension_reduction_name = "umap", 
                            dims_to_plot = c(1, 2), 
                            split_variable = NULL, 
                            plot_axis_lines = TRUE,
                            figure_name_pref = NULL){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  if (missing(split_variable)){
    if (plot_axis_lines == TRUE){
      p <- DimPlot(seuratObject, reduction = dimension_reduction_name, dims = dims_to_plot, label = T, split.by = split_variable, pt.size = 0.5, repel = T, label.size = 10) + 
        theme(
          axis.title.x = element_text(size = 18, face = "bold"), 
          axis.title.y = element_text(size = 18, face = "bold"), 
          axis.ticks.length = unit(.30, "cm"), 
          axis.text = element_text(size = 18, face = "bold"),
          legend.key.size = unit(4,"line"),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.key = element_rect(size = 20),
          legend.spacing = unit(0.05, 'cm'),
          legend.text = element_text(size = 18, face = "bold"),
          legend.spacing.x = unit(0.1, 'cm')) + 
        guides(colour = guide_legend(override.aes = list(size = 8)))
      
      ggsave(plot = p, filename = str_c(temp_dir, toupper(dimension_reduction_name), "_", "all_clusters_dims_", paste(dims_to_plot, collapse = "_"), figure_name_pref, ".png"), width = 14, height = 14, dpi = 300)
    }
    
    else {
      p <- DimPlot(seuratObject, reduction = dimension_reduction_name, dims = dims_to_plot, label = T, split.by = split_variable, pt.size = 0.5, repel = T, label.size = 10) + 
        theme(
          line = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.key.size = unit(4,"line"),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.key = element_rect(size = 20),
          legend.spacing = unit(0.05, 'cm'),
          legend.text = element_text(size = 18, face = "bold"),
          legend.spacing.x = unit(0.1, 'cm')) + 
        guides(colour = guide_legend(override.aes = list(size = 8)))
      
      ggsave(plot = p, filename = str_c(temp_dir, toupper(dimension_reduction_name), "_", "all_clusters_dims_", paste(dims_to_plot, collapse = "_"), figure_name_pref, ".png"), width = 14, height = 14, dpi = 300)
    }
  }
  
  else {
    if (plot_axis_lines == TRUE){
      p <- DimPlot(seuratObject, reduction = dimension_reduction_name, dims = dims_to_plot, label = T, split.by = split_variable, pt.size = 0.5, repel = T, label.size = 10) + 
        theme(
          axis.title.x = element_text(size = 18, face = "bold"), 
          axis.title.y = element_text(size = 18, face = "bold"), 
          axis.ticks.length = unit(.30, "cm"), 
          axis.text = element_text(size = 18, face = "bold"),
          legend.key.size = unit(4,"line"),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.key = element_rect(size = 20),
          legend.spacing = unit(0.05, 'cm'),
          legend.text = element_text(size = 18, face = "bold"),
          legend.spacing.x = unit(0.1, 'cm')) + 
        guides(colour = guide_legend(override.aes = list(size = 8)))
      
      ggsave(plot = p, filename = str_c(temp_dir, toupper(dimension_reduction_name), "_", "all_clusters_dims_", paste(dims_to_plot, collapse = "_"), "_splitBy_", split_variable, figure_name_pref, ".png"), width = 32, height = 18, dpi = 300)
    }
    
    else {
      p <- DimPlot(seuratObject, reduction = dimension_reduction_name, dims = dims_to_plot, label = T, split.by = split_variable, pt.size = 0.5, repel = T, label.size = 10) + 
        theme(
          line = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.key.size = unit(4,"line"),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.key = element_rect(size = 20),
          legend.spacing = unit(0.05, 'cm'),
          legend.text = element_text(size = 18, face = "bold"),
          legend.spacing.x = unit(0.1, 'cm')) + 
        guides(colour = guide_legend(override.aes = list(size = 8)))
      
      ggsave(plot = p, filename = str_c(temp_dir, toupper(dimension_reduction_name), "_", "all_clusters_dims_", paste(dims_to_plot, collapse = "_"), "_splitBy_", split_variable, figure_name_pref, ".png"), width = 32, height = 18, dpi = 300)
    }
  }
}
