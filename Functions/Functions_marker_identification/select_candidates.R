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