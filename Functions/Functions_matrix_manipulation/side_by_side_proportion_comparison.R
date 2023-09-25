# Species bar plot with proportion view per cluster
side_by_side_proportion_comparison <- function(seuratObject, 
                                               store_dir = NULL, 
                                               store_folder = "Proportion_n_count_results",
                                               split_variable_name = "Species", 
                                               figure_name_pref = NULL){
  
  # Manually created color palette to have distinguishable color for the clusters
  mypal = c("#008B45FF", "#075149FF", "#1A9993FF", "#71D0F5FF", "#197EC0FF", "#155F83FF", "#24325FFF", "#631879FF", "#A20056FF", "#FD8CC1FF", 
            "#FB6467FF", "#F05C3BFF", "#BB0021FF", "#800000FF", "#91331FFF", "#C16622FF", "#FFA319FF", "#FED439FF", "#767676FF", "#917C5DFF", 
            "#02401B", "#81A88D", "#354823", "#FAD77B", "#D8B70A", "#A2A475", "#D8A499", "#9986A5", "#CCBA72", "#D9D0D3", "#EAD3BF", 
            "#B6854D", "#F1BB7B", "#FAEFD1", "#74A089", "#CDC08C")
  
  # Creating necessary storing space to store the results
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  if (length(levels(Idents(seuratObject))) <= length(mypal)){
    sps_col = mypal[1:length(levels(Idents(seuratObject)))]
  }
  
  # lets create bar graph for species
  # Let's transform the metadata column into factor
  seuratObject@meta.data[[split_variable_name]] = as.factor(seuratObject@meta.data[[split_variable_name]])
  
  cell_prop_sps = as.data.frame(table(seuratObject@meta.data[[split_variable_name]], Idents(seuratObject)))
  colnames(cell_prop_sps) <- c("species", "cluster", "cell_count")
  
  Species_proportion = cell_prop_sps %>% group_by(species) %>% mutate(rel_prop = cell_count / sum(cell_count))
  
  # Total count of cells from each replicate in each cluster
  p <- ggplot(data = Species_proportion, aes(x = species, y = rel_prop)) + 
    geom_bar(stat = "identity", position = "fill", aes(fill = fct_rev(cluster))) +
    xlab("Species") + ylab("Proportion of cells in each cluster") + 
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "#71D0F5FF"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), # Background of the entire plot
      axis.title = element_text(size = 24, face = "bold", color = "black"),
      axis.ticks.length = unit(.20, "cm"),
      axis.text.x = element_text(size = 24, face = "bold.italic", colour = "black"),
      axis.text.y = element_text(size = 24, face = "bold", colour = "black"),
      legend.key.size = unit(2, "line"), 
      legend.text = element_text(size = 22, face = "bold.italic"),
      legend.title = element_text(size = 22, face = "bold")) + 
    scale_fill_manual(values = mypal) + guides(fill = guide_legend(title = "Clusters"), color = guide_legend(override.aes = list(size = 8)))
  
  ggsave(filename = str_c(temp_dir, figure_name_pref, "Species_n_cluster_proportion.png"), plot = p, width = 12, height = 20, dpi = 300)
}
