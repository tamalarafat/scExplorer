map_factor_loadings <- function(seuratObject,
                                store_dir = NULL,
                                store_folder = "Reduced_dimension_visualization",
                                factor_ids = NULL,
                                factor_to_plot = "inmf",
                                reduction_name = "umap",
                                dims_to_plot = c(1, 2))
  {
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "Factors_loading_on_", toupper(reduction_name)))){
    dir.create(str_c(store_dir, "/", store_folder, "/", "Factors_loading_on_", toupper(reduction_name)), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/", "Factors_loading_on_", toupper(reduction_name), "/")
  
  if (missing(factor_ids)){
    factor_ids = c(1:ncol(seuratObject@reductions[[factor_to_plot]]@cell.embeddings))
  }
  
  for (i in c(1:length(factor_ids))){
    seuratObject$temp_factor = seuratObject@reductions[[factor_to_plot]]@cell.embeddings[, factor_ids[i]]
    
    seuratObject$temp_factor = (seuratObject$temp_factor)*10
    
    p <- FeaturePlot(seuratObject, features = "temp_factor", reduction = reduction_name, dims = dims_to_plot, pt.size = 0.5, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[8])) +
          ggtitle(str_c("GEP", " - " , factor_ids[i])) + labs(color = "Loadings") + 
          theme(
            axis.title = element_text(size = 18, face = "bold"),
            axis.ticks.length = unit(.30, "cm"), 
            axis.text = element_text(size = 18, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.key.size = unit(2,"line"),
            legend.key = element_rect(size = 20),
            legend.text = element_text(size = 18, face = "bold"))
    
    ggsave(filename = str_c(temp_dir, "GEP", "_" , factor_ids[i], "_loading_across_cells_", toupper(reduction_name), ".png"), plot = p, width = 14, height = 14, dpi = 300)
    
  }
}

  
  