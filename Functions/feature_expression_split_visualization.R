# Plot expression of a gene in the reeduced dimension and split the plot by group/replicate/species
feature_expression_split_visualization <- function(seuratObject, 
                                                   store_dir = NULL, 
                                                   store_folder = "Features_plot",
                                                   gene_ID, 
                                                   gene_name = NULL, 
                                                   split_variable = "Species", 
                                                   reduction_name = "umap",
                                                   dims_to_plot = c(1, 2),
                                                   ncol = 2, 
                                                   nrow = 1){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "Genes_expression_split_by_", split_variable))){
    dir.create(str_c(store_dir, "/", store_folder, "/", "Genes_expression_split_by_", split_variable), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/", "Genes_expression_split_by_", split_variable, "/")
  
  
  temp_list = list()
  
  if (missing(gene_name)){
    gene_name = gene_ID
  }
  
  base <- FeaturePlot(seuratObject, features = c(gene_ID), reduction = reduction_name, pt.size = 1, split.by = split_variable, order = T, combine = FALSE)
  
  for (i in c(1:length(base))){
    p <- base[[i]] + theme(axis.text.x = element_text(size = 24, face = "bold"),
                           axis.line.x.bottom = element_line(colour = "darkblue", size = 1, linetype = "solid"), 
                           axis.line.y.left = element_line(colour = "darkblue", size = 1, linetype = "solid"),
                           panel.border = element_blank(),
                           line = element_blank(), 
                           axis.title.x = element_text(size = 24, face = "bold"),
                           axis.title = element_text(size = 24, face = "bold"),
                           axis.ticks.length = unit(.30, "cm"), 
                           axis.text = element_text(size = 24, face = "bold"),
                           legend.key.size = unit(2,"line"),
                           legend.key = element_rect(size = 20),
                           legend.text = element_text(size = 24, face = "bold"),
                           plot.title = element_blank(),
                           legend.title = element_text(size = 24, face = "bold")) +  
      scale_color_gradient(low = "gray", high = RColorBrewer::brewer.pal(9, "Reds")[7], limits = c(0, max(seuratObject@assays$RNA@data[gene_ID, ])),
                           name = gene_name)
    
    temp_list[[i]] <- p
  }
  
  final_fig <- ggarrange(plotlist = temp_list, common.legend = TRUE, ncol = ncol, nrow = nrow, legend = "right")
  
  ggsave(filename = str_c(temp_dir, gene_name, "_expression_across_cells_split_by_", split_variable, "_", toupper(reduction_name), "_dims_", paste(dims_to_plot, collapse = "_"), ".png"), plot = final_fig, 
         width = 32, height = 18, dpi = 300, bg = "white")
  
}
