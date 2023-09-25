# Only for two genes :: plot expression of the genes in a reduced dimension plot
genes_feature_plot <- function(seuratObject,
                               store_dir = NULL,
                               store_folder = "Features_plot",
                               gene_IDs,
                               genes_name = NULL,
                               reduction_name = "umap",
                               dims_to_plot = c(1, 2),
                               figure_name_pref = NULL
                               ){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  
  
  if (length(gene_IDs) == 1){
    p <-  FeaturePlot(seuratObject, features = gene_IDs, reduction = reduction_name, dims = dims_to_plot, pt.size = 1, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[8])) + 
      ggtitle(str_c(gene_IDs, " - " , genes_name)) + 
      theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.ticks.length = unit(.30, "cm"), 
        axis.text = element_text(size = 18, face = "bold"),
        title = element_text(size = 24, face = "bold"),
        legend.key.size = unit(2,"line"),
        legend.key = element_rect(size = 20),
        legend.text = element_text(size = 18, face = "bold"))
    ggsave(filename = str_c(temp_dir, gene_IDs, "_expression_across_cells_", toupper(reduction_name), "_dims_", paste(dims_to_plot, collapse = "_"), figure_name_pref, ".png"), plot = p, width = 14, height = 14, dpi = 300, bg = "transparent")
  }
  
  else {
    p1 <-  FeaturePlot(seuratObject, features = gene_IDs[1], reduction = reduction_name, dims = dims_to_plot, pt.size = 1, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[8])) + 
      ggtitle(str_c(gene_IDs[1], " - " , genes_name[1])) + 
      theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.ticks.length = unit(.30, "cm"), 
        axis.text = element_text(size = 18, face = "bold"),
        title = element_text(size = 24, face = "bold"),
        legend.key.size = unit(2,"line"),
        legend.key = element_rect(size = 20),
        legend.text = element_text(size = 18, face = "bold"))
    
    p2 <-  FeaturePlot(seuratObject, features = gene_IDs[2], reduction = reduction_name, pt.size = 1, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[8])) + 
      ggtitle(str_c(gene_IDs[2], " - " , genes_name[2])) + 
      theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.ticks.length = unit(.30, "cm"), 
        axis.text = element_text(size = 18, face = "bold"),
        title = element_text(size = 24, face = "bold"),
        legend.key.size = unit(2,"line"),
        legend.key = element_rect(size = 20),
        legend.text = element_text(size = 18, face = "bold"))
    
    arranged_fig <- ggarrange(p1, p2, ncol = 2, nrow = 1)
    my_fig <- annotate_figure(arranged_fig)
    ggsave(filename = str_c(temp_dir, gene_IDs[1], "_", gene_IDs[2], "_expression_across_cells_", toupper(reduction_name), "_dims_", paste(dims_to_plot, collapse = "_"), figure_name_pref, ".png"), my_fig,width = 32, height = 18, dpi = 300)
  }
}
