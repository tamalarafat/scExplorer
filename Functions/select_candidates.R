select_candidates <- function(marker_file, find_candidates = 10) {
  
  n_pct = colnames(marker_file)[str_detect(colnames(marker_file), pattern = "pct.2")]
  
  if (!"gene_ID" %in% colnames(cluster_degs)){
    marker_file$gene_ID = rownames(marker_file)
  }
  
  
  if (length(n_pct) == 1){
    # Let's define the candidates selection criteria
    marker_file = marker_file[order(marker_file[["pct.2"]], decreasing = FALSE), ]
    
    if (nrow(marker_file) == 0) {
      print("No TFs")
      } 
    
    else if (nrow(marker_file) <= find_candidates & nrow(marker_file) != 0){
      marker_file = marker_file
    } 
    
    else {
      if (nrow(marker_file[marker_file$avg_log2FC >= 1,]) < find_candidates & nrow(marker_file[marker_file$avg_log2FC >= 1,]) > 0){
        subset1 = marker_file[marker_file$avg_log2FC >= 1,]
        subset2 = marker_file[marker_file$avg_log2FC < 1,]
        subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
        
        if (nrow(subset2) <= (find_candidates - nrow(subset1))) {
          marker_file = rbind.data.frame(subset1, subset2)
        } else {
          subset2 = head(subset2, (find_candidates - nrow(subset1)))
          marker_file = rbind.data.frame(subset1, subset2)
        }
      }
      else if (nrow(marker_file[marker_file$avg_log2FC >= 1,]) == 0){
        marker_file = marker_file[order(marker_file$avg_log2FC, decreasing = TRUE), ]
        marker_file = marker_file[c(1:find_candidates), ]
      }
      else {
        marker_file = marker_file[marker_file$avg_log2FC >= 1, ]
        marker_file = marker_file[c(1:find_candidates), ]
      }
    }
    
  }
  
  else if (length(n_pct) == 2){
    marker_file = marker_file[order(marker_file[, n_pct[1]], marker_file[, n_pct[2]]), ]
    rownames(marker_file) = marker_file$gene_ID
    
    avg_FC_names = colnames(marker_file)[str_detect(colnames(marker_file), pattern = "avg_log2FC")]
    
    if (nrow(marker_file) == 0) {
      print("No TFs")} 
    
    else if (nrow(marker_file) <= find_candidates & nrow(marker_file) != 0){
      marker_file = marker_file
    } 
    
    else {
      if (nrow(marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]) < find_candidates & 
          nrow(marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]) > 0){
        subset1 = marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]
        subset2 = marker_file[marker_file[[avg_FC_names[1]]] < 1 & marker_file[[avg_FC_names[2]]] < 1, ]
        subset2 = subset2[order(subset2[, avg_FC_names[1]], subset2[, avg_FC_names[2]], decreasing = TRUE), ]
        
        if (nrow(subset2) <= (find_candidates - nrow(subset1))) {
          marker_file = rbind.data.frame(subset1, subset2)
        } 
        
        else {
          subset2 = head(subset2, (find_candidates - nrow(subset1)))
          marker_file = rbind.data.frame(subset1, subset2)
        }
      }
      else if (nrow(marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]) == 0){
        marker_file = marker_file[order(marker_file[, avg_FC_names[1]], marker_file[, avg_FC_names[2]], decreasing = TRUE), ]
        marker_file = marker_file[c(1:find_candidates), ]
      }
      
      else {
        marker_file = marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]
        marker_file = marker_file[c(1:find_candidates), ]
      }
    }
  }
  
  marker_file$source = "DEGs"
  
  return(marker_file)
}


# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("~/Documents/2023/PhD_projects_Yasir/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

conserved_marker = loadRData("/home/yasir/Documents/Projects_Yasir/comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs/Differentially_expressed_genes/Conserved_DEGs_grouped_by_Species/Conserved_markers_DEtest_wilcox/Cluster_15_positive_conserved_marker.RData")

cluster_degs = loadRData("/home/yasir/Documents/Projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox/Cluster_15_DEGs.RData")

scm = specific_conserved_marker_finder(DEG_file = conserved_marker, pct_diff = 0.3, include_pct_diff = TRUE, store_outputs = FALSE)

sm = specific_marker_finder(DEG_file_dir = cluster_degs, max_pct2_detection = 0.1, pct_diff = 0.3, store_outputs = FALSE)

ccm = select_candidates(marker_file = scm, find_candidates = 10)

csm = select_candidates(marker_file = sm, find_candidates = 10)

# Path to the DEG files - differentially genes usign the "Findmarkers" function (Seurat) per cluster
DEG_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

candidate_markers_DEGs(DEG_file = DEG_dir, file_name_pattern = "Cluster_", max_pct2_detection = 0.1, pct_diff = 0.3, include_pct_diff = TRUE, find_candidates = 10, store_outputs = TRUE)

# Load the TFs list 
ch_TFs = read.delim("/home/ytamal2/Documents/2023/PhD_projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/Analysis_objects/Annotation_files/TFs_list/cardamine_TF_IDs.txt", header = FALSE)[, 1]

ml = candidate_markers_DEGs(DEG_file = DEG_dir, file_name_pattern = "Cluster_", max_pct2_detection = 0.1, pct_diff = 0.3, include_pct_diff = TRUE, find_candidates = 10, store_outputs = FALSE, specify_gene_ID = ch_TFs, incorporate_column_name = "TF", combine_categories = TRUE)

ml_15 = candidate_markers_DEGs(DEG_file = DEG_dir, file_name_pattern = "Cluster_", max_pct2_detection = 0.1, pct_diff = 0.3, include_pct_diff = TRUE, find_candidates = 10, store_outputs = FALSE, specify_gene_ID = ch_TFs, incorporate_column_name = "TF", combine_categories = TRUE, cluster_ID = 15)
