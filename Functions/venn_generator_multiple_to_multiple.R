
# Load and store the csv files
venn_generator_multiple_to_multiple <- function(First_marker_list, Second_marker_list){
  
  # Let's create a directory to store the figures if the directory does not exit
  if (!dir.exists("Venn_Figures_markers_comparison")){
    dir.create("Venn_Figures_markers_comparison", showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # Lets get the names of the DEG files in the list so we can iterate over the list items
  temp_file_names_1 = names(First_marker_list)
  temp_file_names_2 = names(Second_marker_list)
  
  for (i in c(1:length(temp_file_names_1))){
    
    if (!dir.exists(str_c("Venn_Figures_markers_comparison/", temp_file_names_1[i]))){
      dir.create(str_c("Venn_Figures_markers_comparison/", temp_file_names_1[i]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Get the DEG file stored as item in the list
    perm_DEG = First_marker_list[[temp_file_names_1[i]]] # Let's say it's the permanent item which we want to compare with all other items in the second list
    
    # Now let's iterate over the second lsit
    for (j in c(1:length(temp_file_names_2))){
      temp_DEG = Second_marker_list[[temp_file_names_2[j]]] # Let's say it's the temporary item which will change in each iteration
      
      # Let's create a list with the two items (i, j)
      vennInput = list(L1 = rownames(perm_DEG), L2 = rownames(temp_DEG))
      names(vennInput) = c(str_c("List_1_", temp_file_names_1[i]), str_c("List_2_", temp_file_names_2[j]))
      
      # Generate the venn diagram for the list components
      p <- ggvenn(vennInput, columns = c(str_c("List_1_", temp_file_names_1[i]), str_c("List_2_", temp_file_names_2[j])), fill_color = c("#075149FF", "#FFA319FF"), 
                  fill_alpha = 0.2, text_size = 8, set_name_size = 8, stroke_alpha = 0.2, stroke_size = 0.2)
      
      ggsave(filename = str_c("Venn_Figures_markers_comparison/", temp_file_names_1[i], "/List_1_", temp_file_names_1[i], "_", "List_2_", temp_file_names_2[j], ".png"), plot = p, width = 18, height = 18, dpi = 300, bg = "white")
    }
  }
}
