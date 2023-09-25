# Function to generate a gene's count and proportion per cluster
feature_count_n_proportion <- function(seuratObject, 
                                       store_dir = NULL, 
                                       store_folder = "Proportion_n_count_results",
                                       gene_ID,
                                       gene_name,
                                       figure_name = NULL){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # Lets create a empty vector to store the count information
  temp = c()
  
  clusterLabels = str_sort(levels(Idents(seuratObject)), numeric = TRUE)
  
  for (i in c(1:length(clusterLabels))){
    cellNames = WhichCells(seuratObject, idents = clusterLabels[i])
    cellCount = sum(GetAssayData(seuratObject, assay = "RNA", slot = "counts")[gene_ID, cellNames] != 0)
    temp = c(temp, cellCount)
  }
  
  cell_prop = as.data.frame(table(Idents(seuratObject)))
  colnames(cell_prop) <- c("cluster", "total_cell")
  
  cell_prop = cell_prop %>% mutate(detected = temp) %>% 
    mutate(non_detected = total_cell - detected) %>% 
    mutate(frac_detected = round((detected / sum(detected)) * 100, 1)) %>% 
    mutate(frac_label = str_c(frac_detected, "%")) %>% 
    mutate(count_label = str_c(detected, frac_label, sep = ",\n"))
  
  
  long_cell_prop = cell_prop %>% dplyr::select(-c(frac_detected, frac_label, count_label)) %>% 
    gather("detection", "cell_count", detected, non_detected) %>% 
    group_by(cluster) %>% 
    mutate(cluster_prop = round((cell_count/total_cell) * 100, 2)) %>% 
    mutate_at(vars(detection), as.factor)
  
  # frac_detected = In each cluster, the proportion of cells in which expression of the given gene is detected (cells with the genes expression in a cluster / total cells with that genes expression)
  # cluster_prop = ration of detected and non detected cells in each cluster
  
  p1 <- ggplot(data = long_cell_prop, aes(x = fct_rev(cluster), y = cell_count)) + 
    geom_bar(stat = "identity", position = "stack", aes(fill = detection), width = 0.6) + 
    geom_text(data = cell_prop, aes(x = cluster, y = (total_cell + 500), label = count_label),  size = 6, fontface = "bold") + 
    coord_flip() + 
    xlab("Cell clusters") + ylab("Number of cells (cells in which expression was detected, proportion of detected cells)") + 
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "#71D0F5FF"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), # Background of the entire plot
      axis.title = element_text(size = 22, face = "bold", color = "black"),
      axis.ticks.length = unit(.20, "cm"),
      axis.text = element_text(size = 22, face = "bold", colour = "black"),
      legend.key.size = unit(2, "line"), 
      legend.text = element_text(size = 22, face = "bold"),
      legend.title = element_text(size = 22, face = "bold")) + 
    scale_fill_manual(values = c("#00A08A", "#F98400"), breaks = c("non_detected", "detected"), labels = c("not detected", "detected")) + 
    guides(fill = guide_legend(title = "Expression"))
  
  p2 <- ggplot(data = long_cell_prop, aes(x = fct_rev(cluster), y = cluster_prop)) + 
    geom_bar(stat = "identity", position = "fill", aes(fill = detection), width = 0.6) + 
    coord_flip() + 
    xlab("Cell clusters") + ylab("Fraction of cells with detected expression") + 
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "#71D0F5FF"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), # Background of the entire plot
      axis.title = element_text(size = 22, face = "bold", color = "black"),
      axis.ticks.length = unit(.20, "cm"),
      axis.text = element_text(size = 22, face = "bold", colour = "black"),
      legend.key.size = unit(2, "line"), 
      legend.text = element_text(size = 22, face = "bold"),
      legend.title = element_text(size = 22, face = "bold")) + 
    scale_fill_manual(values = c("#00A08A", "#F98400"), breaks = c("non_detected", "detected"), labels = c("not detected", "detected")) + 
    guides(fill = guide_legend(title = "Expression"))
  
  arranged_fig <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
  ggsave(filename = str_c(temp_dir, gene_name, "_count_n_proportion", figure_name, ".png"), plot = arranged_fig, width = 32, height = 22, dpi = 300, bg = "white")
}
