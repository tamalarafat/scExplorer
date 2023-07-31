# Load and store the csv files
venn_generator_multiple_to_one <- function(first_marker_list, second_marker_list){
  
  # Let's create a directory to store the figures if the directory does not exit
  if (!dir.exists("Venn_diagram_of_overlapping_genes")){
    dir.create("Venn_diagram_of_overlapping_genes", showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c("Venn_diagram_of_overlapping_genes", "/")
  
  # Lets get the names of the DEG files in the list so we can iterate over the list items
  if (!is.list(first_marker_list)) {
    first_marker_list = list(first_marker_list)
    names(first_marker_list)[1] = "marker_set_1"
  }
  
  if (!is.list(second_marker_list)) {
    second_marker_list = list(second_marker_list)
    names(second_marker_list)[1] = "marker_set_2"
  }
  
  temp_file_names_1 = names(first_marker_list)
  temp_file_names_2 = names(second_marker_list)
  
  for (i in c(1:length(temp_file_names_1))){
    # Get the DEG file stored as item in the list
    first_list = first_marker_list[[temp_file_names_1[i]]] # Let's say it's the permanent item which we want to compare with all other items in the second list
    
    # Now let's iterate over the second lsit
    for (j in c(1:length(temp_file_names_2))){
      second_list = second_marker_list[[temp_file_names_2[j]]] # Let's say it's the temporary item which will change in each iteration
      
      # Let's create a list with the two items (i, j)
      venn_input = list(L1 = first_list, L2 = second_list)
      names(venn_input) = c(temp_file_names_1[i], temp_file_names_2[j])
      
      # Generate the venn diagram for the list components
      p <- ggvenn(venn_input, columns = c(temp_file_names_1[i], temp_file_names_2[j]), fill_color = c("#075149FF", "#FFA319FF"), 
                  fill_alpha = 0.2, text_size = 8, set_name_size = 8, stroke_alpha = 0.2, stroke_size = 0.2)
      
      ggsave(filename = str_c(temp_dir, temp_file_names_1[i], "_", temp_file_names_2[j], ".png"), plot = p, width = 18, height = 18, dpi = 300, bg = "white")
    }
  }
}