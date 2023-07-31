
# Calculate specificity score of the GEPs / latent factors per group (Genotypes or Species)
calculate_GEP_specificity_for_groups <- function(LigerObject){
  # What are the datasets
  dfs <- levels(LigerObject@cell.data$dataset)
  
  # Lets remove digits from strings, so we have only the species/genotypes initials
  orgs <- unique(gsub("[0-9]+", "", dfs))
  
  # if we have more than two organisms, we need to define that many number of empty variables
  org1 = dfs[str_detect(dfs, pattern = orgs[1])] # output is a boolean for the whole string
  org2 = dfs[grep(pattern = orgs[2], dfs)] # output is a boolean for the whole string
  
  # Let's create a directory to store the figures
  
  for (i in c(1:length(org1))){
    for (j in c(1:length(org2))){
      DS <- calcDatasetSpecificity(LigerObject, dataset1 = org1[i], dataset2 = org2[j], do.plot = F)
      
      if (!exists("DF")){
        DF <- data.frame(Specificity = DS[[3]], Factors = c(1:length(DS[[3]])))
        DF <- DF[ ,c(2,1)]
        DF$Factors <- factor(DF$Factors, levels = seq(1, length(DF$Factors)))
      }
      else {
        temp_df <- data.frame(Specificity = DS[[3]], Factors = c(1:length(DS[[3]])))
        DF <- cbind.data.frame(DF, temp_df$Specificity)
      }
    }
  }
  DF$MeanSpecificity = rowMeans(DF[, -c(1)])
  DF = DF[ ,c("Factors", "MeanSpecificity")]
  return(DF)
} 
