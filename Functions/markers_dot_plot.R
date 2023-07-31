# Generate dotplot - expression of marker genes across the cell clusters
markers_dot_plot <- function(seuratObject, 
                             store_dir = NULL, 
                             store_folder = "Markers_dotplot",
                             marker_file, 
                             genes_ID_column, 
                             genes_name_column, 
                             split_variable_name = NULL, 
                             group_clusters = FALSE, 
                             figure_name_pref = NULL){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # To plot the genes with the row names, here I'm assigning the gene_ID column name as the row name of the marker file
  rownames(marker_file) <- marker_file[ , genes_ID_column]
  
  if (missing(split_variable_name)){
    p <- DotPlot(seuratObject, features = rownames(marker_file), scale = T, assay = "RNA", cols = RColorBrewer::brewer.pal(9, name = "Greens")[c(2,9)], col.min = 0, dot.scale = 25, dot.min = .05, cluster.idents = group_clusters) +
      guides(size = guide_legend(override.aes = list(color= "orange", alpha = 1))) +
      xlab("Genes") +
      ylab("Clusters") +
      scale_x_discrete(labels = as.character(marker_file[, genes_name_column])) + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 90),
            axis.title.x = element_text(size = 24, face = "bold"),
            axis.title.y = element_text(size = 24, face = "bold"),
            axis.ticks.length = unit(.20, "cm"),
            axis.text = element_text(size = 24, face = "bold", colour = "black"),
            legend.key.size = unit(2,"line"),
            legend.text = element_text(size = 24, face = "bold"),
            legend.spacing = unit(2.0, 'cm'), 
            legend.position = "top", 
            legend.key.width = unit(2.0, 'cm'),
            legend.title = element_text(size = 24, colour = "black", face = "bold"))
    
    ggsave(filename = str_c(temp_dir, figure_name_pref, "known_cell_type_markers_expression.png"), plot = p, width = 32, height = 22, dpi = 300)
  }
  
  else {
    p <- DotPlot(seuratObject, features = rownames(marker_file), dot.scale = 14, cluster.idents = group_clusters, dot.min = .05, split.by = split_variable_name, 
                 cols = c(RColorBrewer::brewer.pal(9, "Greens")[7] , RColorBrewer::brewer.pal(9, "Reds")[7])) + RotatedAxis() +
      scale_x_discrete(labels = as.character(marker_file[, genes_name_column])) +
      xlab("Genes") +
      ylab("Split clusters identity") +
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 90),
            axis.title.x = element_text(size = 24, face = "bold"),
            axis.title.y = element_text(size = 24, face = "bold"),
            axis.ticks.length = unit(.20, "cm"),
            axis.text = element_text(size = 24, face = "bold", colour = "black"),
            legend.key.size = unit(2,"line"),
            legend.text = element_text(size = 24, face = "bold"),
            legend.spacing = unit(2.0, 'cm'), 
            legend.position = "top", 
            legend.key.width = unit(2.0, 'cm'),
            legend.title = element_text(size = 24, colour = "black", face = "bold"))
    
    ggsave(filename = str_c(temp_dir, figure_name_pref, "known_cell_type_markers_expression_split_by_", split_variable_name, ".png"), plot = p, width = 32, height = 36, dpi = 300)

  }
}
