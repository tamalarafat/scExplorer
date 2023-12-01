was_it_installed <- function(package_name, retreive_info = "Version"){
  
  # package_name - name of the package you're interested in
  
  # Get information about installed packages
  installed_packages <- installed.packages()
  
  # Check if the package is installed
  if (package_name %in% installed_packages[, "Package"]) {
    # Print the version of the package
    print(paste("Yes, installed. Version of", package_name, ":", installed_packages[package_name, retreive_info]))
  } 
  
  else {
    print(paste(package_name, "is not installed."))
  }
  
}
