# Function to write text files
fileGenerator <- function(markerList, 
                          fileName) {
  
  # initialize the connection to the file name (specified by the fileName) 
  fileWriter <- file(fileName)
  
  # writes the contents of the markerList to the file using the file connection fileWriter
  writeLines(markerList, fileWriter)
  
  # Close the file connection
  close(fileWriter)
}
