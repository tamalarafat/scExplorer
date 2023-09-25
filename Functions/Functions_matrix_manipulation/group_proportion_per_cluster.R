# The function name has been changed from species_proportion_per_cluster to group_proportion_per_cluster - 24/03/2023
# As we will group the replicates based on the species, we need both species information and replicates information
group_proportion_per_cluster <- function(seuratObject, 
                                    store_dir = NULL, 
                                    store_folder = "Proportion_n_count_results",
                                    replicate_metadata_name = "Datasets",
                                    split_metadata_name = "Species",
                                    rep_prop = FALSE){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # Let's transform the metadata column into factor
  seuratObject@meta.data[[replicate_metadata_name]] = as.factor(seuratObject@meta.data[[replicate_metadata_name]])
  
  # Let's transform the metadata column into factor - changed 13/03/2023
  if (class(seuratObject@meta.data[[split_metadata_name]]) != "factor"){
    seuratObject@meta.data[[split_metadata_name]] = as.factor(seuratObject@meta.data[[split_metadata_name]])
  }
  
  # Let's create a data table with cell count and cluster idents
  cell_prop_sps = as.data.frame(table(seuratObject@meta.data[[split_metadata_name]], Idents(seuratObject)))
  colnames(cell_prop_sps) <- c("species", "cluster", "cell_count")
  
  temp_sps = levels(seuratObject@meta.data[[split_metadata_name]])
  
  temp_sps_list = list()
  
  for (i in c(1:length(temp_sps))){
    temp_df = cell_prop_sps[cell_prop_sps$species == temp_sps[i], ]
    temp_df$sps_prop = round((temp_df$cell_count/sum(temp_df$cell_count)) * 100, 2)
    temp_df$prop_label = paste(as.character(temp_df$sps_prop), "%", sep = "")
    temp_sps_list[[i]] = temp_df
  }
  
  Species_proportion = do.call(rbind.data.frame, temp_sps_list)
  Species_proportion$cluster <- factor(Species_proportion$cluster, levels = str_sort(levels(Species_proportion$cluster), numeric = TRUE))
  rownames(Species_proportion) <- NULL
  
  
  # set the color base
  base_col = "#f2edee" # The base color
  grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A") # Manually picked colors for each species class
  sps_col = grp_col[c(1:length(temp_sps))]
  
  if (rep_prop == FALSE){
    p <- ggplot(data = Species_proportion, aes(x = cluster, y = sps_prop)) + 
      geom_bar(stat = "identity", position = "fill", aes(fill = species)) +
      xlab("Clusters") + ylab("Proportion of cells") + 
      coord_flip() + 
      theme(
        panel.border = element_blank(),
        axis.line = element_line(colour = "#71D0F5FF"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), # Background of the entire plot
        axis.title = element_text(size = 24, face = "bold", color = "black"),
        axis.ticks.length = unit(.20, "cm"),
        axis.text = element_text(size = 24, face = "bold", colour = "black"),
        legend.key.size = unit(2, "line"), 
        legend.text = element_text(size = 22, face = "bold.italic"),
        legend.title = element_text(size = 22, face = "bold")) + 
      scale_fill_manual(values = sps_col) + guides(fill = guide_legend(title = "Species"), color = guide_legend(override.aes = list(size = 8)))
    
    ggsave(filename = str_c(temp_dir, "Proportion_of_cells_of_", split_metadata_name, "_per_cluster.png"), plot = p, width = 22, height = 28, dpi = 300)
  }
  
  else {
    
    # Let's create a data table with cell count and cluster idents
    cell_prop_rep = as.data.frame(table(seuratObject@meta.data[[replicate_metadata_name]], Idents(seuratObject)))
    colnames(cell_prop_rep) <- c("replicates", "cluster", "cell_count")
    
    temp_reps = levels(seuratObject@meta.data[[replicate_metadata_name]])
    
    temp_reps_list = list()
    
    for (i in c(1:length(temp_reps))){
      temp_df = cell_prop_rep[cell_prop_rep$replicates == temp_reps[i], ]
      temp_df$rep_prop = round((temp_df$cell_count/sum(temp_df$cell_count)) * 100, 2)
      temp_df$prop_label = paste(as.character(temp_df$rep_prop), "%", sep = "")
      temp_reps_list[[i]] = temp_df
      
    }
    
    Replicates_proportion = do.call(rbind.data.frame, temp_reps_list)
    Replicates_proportion$cluster <- factor(Replicates_proportion$cluster, levels = str_sort(levels(Replicates_proportion$cluster), numeric = TRUE))
    rownames(Replicates_proportion) <- NULL
    
    # As the replicates proportion is in higher magnitude (multiplied by 100), we need to rescale the data
    rep_df = Replicates_proportion %>% group_by(cluster) %>% mutate(rel_prop = rep_prop / sum(rep_prop)) %>% mutate(sps_group = gsub("[^A-Z]", "", replicates)) %>% mutate(sps_group = as.factor(sps_group))
    
    p <- ggplot(data = Species_proportion, aes(x = cluster, y = sps_prop)) + 
      geom_bar(stat = "identity", position = "fill", aes(fill = species)) +
      xlab("Clusters") + ylab("Proportion of cells") + 
      coord_flip() + 
      theme(
        panel.border = element_blank(),
        axis.line = element_line(colour = "#71D0F5FF"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), # Background of the entire plot
        axis.title = element_text(size = 24, face = "bold", color = "black"),
        axis.ticks.length = unit(.20, "cm"),
        axis.text = element_text(size = 24, face = "bold", colour = "black"),
        legend.key.size = unit(2, "line"), 
        legend.text = element_text(size = 22, face = "bold.italic"),
        legend.title = element_text(size = 22, face = "bold")) + 
      scale_fill_manual(values = sps_col) + guides(fill = guide_legend(title = "Species"), color = guide_legend(override.aes = list(size = 8)))
    
    p <- p + geom_point(data = rep_df, aes(x = cluster, y = rel_prop, color = sps_group, shape = sps_group), size = 8, stroke = 4) + scale_shape_manual(values = c(21, 22)) + 
      scale_color_manual(values = c("#E7B800", "#E7B800")) + guides(shape = guide_legend(title = "Species"), color = guide_legend(title = "Species"))
    
    ggsave(filename = str_c(temp_dir, "Proportion_of_cells_of_", split_metadata_name, "_with_", replicate_metadata_name, "_proportion_per_cluster.png"), plot = p, width = 22, height = 28, dpi = 300)
  }
}
