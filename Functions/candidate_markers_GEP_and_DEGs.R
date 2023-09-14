candidate_markers_GEP_and_DEGs <- function(seurat_object, 
                                           GEP_IDs = NULL, 
                                           cluster_ID,
                                           DEG_file,
                                           DE_test = "wilcox", 
                                           reduction_name = "inmf",
                                           find_candidates = 10,
                                           max_pct2_detection,
                                           pct_diff,
                                           include_pct_diff,
                                           specify_gene_ID,
                                           incorporate_column_name = NULL,
                                           store_outputs = TRUE,
                                           store_dir = NULL, 
                                           store_folder = "Cluster_candidate_biomarkers") {
  
  
  # Creating necessary storing space to store the results
  
  if (store_outputs == TRUE) {
    # Output directory
    if (missing(store_dir)){
      store_dir = getwd()
    }
    
    # Output folder
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Assigning the output directory path to a variable
    temp_dir = str_c(store_dir, "/", store_folder, "/")
  }
  
  if (missing(GEP_IDs)) {
    
    # If GEP ID or IDs are not provided, use the characteristic GEP of the cell cluster to identify candidate markers
    defining_GEP = characteristic_GEP_of_cells(seurat_object = seurat_object, target_ID = cluster_ID, store_output = FALSE)
    
    # Get the ID of the characteristic GEP
    GEP_IDs = as.character(defining_GEP[which.max(defining_GEP$cell_count), "GEP"])
    
  }
  
  # Lets get the 
  gep_candidates = candidate_markers_GEP(seurat_object = seurat_object, 
                                         GEP_IDs = GEP_IDs, 
                                         cluster_ID = cluster_ID, 
                                         find_candidates = find_candidates, 
                                         combine_GEPs = TRUE, 
                                         store_outputs = FALSE)
  
  
  if (!is.null(incorporate_column_name)) {
    gep_candidates[[incorporate_column_name]] = ifelse(gep_candidates$gene_ID %in% specify_gene_ID, "Yes", "No")
  }
  

  
  gep_candidates$source = ifelse(gep_candidates$gene_ID %in% rownames(DEG_file), str_c(gep_candidates$source, "and DEGs", sep = " "), gep_candidates$source)
  
  DEG_file = DEG_file[!rownames(DEG_file) %in% gep_candidates$gene_ID, ]
  
  deg_candidates = candidate_markers_DEGs(DEG_file = DEG_file, 
                                          max_pct2_detection = max_pct2_detection, 
                                          pct_diff = pct_diff, 
                                          include_pct_diff = include_pct_diff, 
                                          find_candidates = find_candidates, 
                                          cluster_ID = cluster_ID, 
                                          specify_gene_ID = specify_gene_ID, 
                                          incorporate_column_name = incorporate_column_name, 
                                          combine_categories = TRUE, 
                                          store_outputs = FALSE)
  
  # Check if the deg table is a typical table for a single species/dataset or a conserved deg table of two species or conditions
  n_pct = colnames(deg_candidates)[str_detect(colnames(deg_candidates), pattern = "pct.2")]
  
  if (length(n_pct) == 1) {
    deg_candidates = deg_candidates
  }
  
  else if (length(n_pct) == 2) {
    pct1_names = colnames(deg_candidates)[str_detect(colnames(deg_candidates), pattern = "pct.1")]
    
    deg_candidates$pct.1 = apply(deg_candidates[, pct1_names], 1, mean)
    
    pct2_names = colnames(deg_candidates)[str_detect(colnames(deg_candidates), pattern = "pct.2")]
    
    deg_candidates$pct.2 = apply(deg_candidates[, pct2_names], 1, mean)
    
    avg_FC_names = colnames(markers_description)[str_detect(colnames(markers_description), pattern = "avg_log2FC")]
    
    deg_candidates$avg_log2FC = apply(deg_candidates[, avg_FC_names], 1, mean)
  }
  
  
  shared_colnames = colnames(gep_candidates)[colnames(gep_candidates) %in% colnames(deg_candidates)]
  
  gep_candidates = gep_candidates[, shared_colnames]
  
  deg_candidates = deg_candidates[, shared_colnames]
  
  cluster_candidates = rbind.data.frame(gep_candidates, deg_candidates)
  
  # Order the data.frame according to the pct.2 criteria where lowest pct corresponds to the top row / order
  cluster_candidates = cluster_candidates[order(cluster_candidates$pct.2, decreasing = FALSE), ]
  
  # Add rank information - ranks are based on the coefficient. Largest value top rank, lowest value bottom rank
  cluster_candidates$rank = c(1:nrow(cluster_candidates))
  
  if (store_outputs == TRUE) {
    # Write csv files for each list item
    write.csv(cluster_candidates, file = str_c(temp_dir, "cluster_", cluster_ID, "_candidate_biomarkers.csv"), row.names = FALSE)
  }
  
  else {
    return(cluster_candidates)
  }
}

