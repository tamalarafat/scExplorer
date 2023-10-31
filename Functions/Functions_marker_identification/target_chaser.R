# Here, add one argument include_pct_diff - this will prompt the user if they want to include the pct difference to select specific markers
target_chaser <- function(DEG_file, 
                          store_dir, 
                          store_folder = "Potential_target_markers", 
                          store_outputs = FALSE, 
                          marker_file_name = "Target_genes.csv"){
  
  # Create necessary storing space to store the results
  
  # Check if the user want to store the outputs. If Yes, create directorries on the desider location
  if (store_outputs == TRUE) {
    
    # If no store directory location in provided, the output folders will be created on the current working directory 
    if (missing(store_dir)){
      store_dir = getwd()
    }
    
    # Provided folder name in which the outputs will be saved
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Variable containing the directory path
    temp_dir = str_c(store_dir, "/", store_folder, "/")
  }
  
  # Get the colnames with "pct" information
  pct1_name = colnames(DEG_file)[str_detect(colnames(DEG_file), pattern = "pct.1")]
  
  pct2_name = colnames(DEG_file)[str_detect(colnames(DEG_file), pattern = "pct.2")]
  
  # Filtering criteria
  potential_markers = DEG_file[(DEG_file[[pct1_name]] - DEG_file[[pct2_name]]) >= DEG_file[[pct2_name]], ]
  
  potential_markers = potential_markers[potential_markers$p_val < 0.01, ]
  
  # checks if the user want to save the markers file in the directory or not
  if (store_outputs == TRUE){
    write.csv(potential_markers, file = str_c(temp_dir, marker_file_name), row.names = TRUE)
  } 
  else {
    potential_markers
  }
  
  }

