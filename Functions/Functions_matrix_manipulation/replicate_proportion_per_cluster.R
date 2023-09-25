# As we will group the replicates based on the species, we need both species information and replicates information
replicate_proportion_per_cluster <- function(seuratObject, 
                                             store_dir = NULL, 
                                             store_folder = "Proportion_n_count_results",
                                             replicate_metadata_name = "Replicates",
                                             split_metadata_name = NULL
){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # Set the assay to RNA
  DefaultAssay(seuratObject) <- "RNA"
  
  # Let's transform the metadata column into factor
  if (class(seuratObject@meta.data[[replicate_metadata_name]]) != "factor"){
    seuratObject@meta.data[[replicate_metadata_name]] = as.factor(seuratObject@meta.data[[replicate_metadata_name]])
  }
  
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
  
  
  # set the color base
  base_col = "#f2edee" # The base color
  grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A") # Manually picked colors for each species class
  
  if (missing(split_metadata_name)){
    rep_col = colorRampPalette(c(base_col, grp_col[2]))(length(temp_reps) + 1)[-1]
    
    p <- ggplot(data = Replicates_proportion, aes(x = cluster, y = rep_prop, fill = replicates)) + 
      geom_bar(stat = "identity", position = "fill") +
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
      scale_fill_manual(values = rep_col) + guides(fill = guide_legend(title = "Replicates"), color = guide_legend(override.aes = list(size = 8)))
    
    ggsave(filename = str_c(temp_dir, "Proportion_of_cells_of_", replicate_metadata_name, "_per_cluster.png"), plot = p, width = 22, height = 28, dpi = 300)
  }
  
  else {
    # Set the color of species if defined while calling the function
    
    # Let's transform the metadata column into factor - changed 23/03/2023 & changed 24/04/2023
    if (class(seuratObject@meta.data[[split_metadata_name]]) != "factor"){
      seuratObject@meta.data[[split_metadata_name]] = as.factor(seuratObject@meta.data[[split_metadata_name]])
    }
    
    # get the species levels
    temp_sps = levels(seuratObject@meta.data[[split_metadata_name]])
    
    sps_col = grp_col[c(1:length(temp_sps))]
    
    temp_col_list = list()
    
    # Lets create color gradient for the replicates for each group of species
    for (items in c(1:length(temp_sps))){
      temp_col = colorRampPalette(c(base_col, sps_col[items]))(length(grep(temp_sps[items], unique(str_c(seuratObject@meta.data[[replicate_metadata_name]], seuratObject@meta.data[[split_metadata_name]])))) + 1)[-1]
      temp_col_list[[items]] = temp_col
    }
    
    Replicates_proportion$groups = gsub(pattern = "[^[:alpha:]]", replacement = "", x = Replicates_proportion$replicates)
    
    temp_groups = unique(Replicates_proportion$groups)
    
    for (i in c(1:length(temp_groups))){
      temp_indc = which(lengths(temp_col_list) == sum(gsub(pattern = "[^[:alpha:]]", replacement = "", x = unique(Replicates_proportion$replicates)) == temp_groups[i]))
      
      reps_per_group = unique(Replicates_proportion[Replicates_proportion$groups == temp_groups[i], ]$replicates)
      
      for (j in c(1:length(reps_per_group))){
        Replicates_proportion[Replicates_proportion$replicates == reps_per_group[j], ]$groups = temp_col_list[[temp_indc]][j]
      }
      temp_col_list = temp_col_list[-c(temp_indc)]
    }
    
    p <- ggplot(data = Replicates_proportion, aes(x = cluster, y = rep_prop, fill = replicates)) + 
      geom_bar(stat = "identity", position = "fill") +
      xlab("Clusters") + ylab("Proportion of cells") + 
      coord_flip() + 
      theme(
        panel.border = element_blank(),
        axis.line = element_line(colour = "#71D0F5FF"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), # Background of the entire plot
        axis.title = element_text(size = 28, face = "bold", color = "black"),
        axis.ticks.length = unit(.20, "cm"),
        axis.text = element_text(size = 28, face = "bold", colour = "black"),
        legend.key.size = unit(2, "line"), 
        legend.text = element_text(size = 28, face = "bold.italic"),
        legend.title = element_text(size = 28, face = "bold")) + 
      scale_fill_manual(values = levels(as.factor(Replicates_proportion$groups))) + 
      guides(fill = guide_legend(title = replicate_metadata_name), color = guide_legend(override.aes = list(size = 8)))
    
    ggsave(filename = str_c(temp_dir, "Proportion_of_cells_of_", replicate_metadata_name, "_per_cluster_by_", split_metadata_name, ".png"), plot = p, width = 22, height = 28, dpi = 300)
  }
} 
