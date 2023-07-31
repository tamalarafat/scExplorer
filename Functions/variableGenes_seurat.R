source("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions/seurat_FindVariableFeatures.default.R")

# Variable gene selection method - How to choose top variable features
variableGenes_seurat <- function(dataObject, 
                                selection.method = "vst", 
                                loess.span = 0.3,
                                clip.max = 'auto',
                                mean.function = FastExpMean, 
                                dispersion.function = FastLogVMR,
                                num.bin = 20,
                                binning.method = "equal_width",
                                mean.cutoff = c(0.1, 8),
                                dispersion.cutoff = c(1, Inf),
                                nfeatures = 3000,
                                verbose = TRUE
){
  # vst :: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). 
  # Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). 
  # Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
  
  # loess.span :: Loess span parameter used when fitting the variance-mean relationship
  
  # clip.max :: After standardization values larger than clip.max will be set to clip.max; default is 'auto' which sets this value to the square root of the number of cells
  
  # mean.function :: Function to compute x-axis value (average expression). Default is to take the mean of the detected (i.e. non-zero) values
  
  # dispersion.function :: Function to compute y-axis value (dispersion). Default is to take the standard deviation of all values
  
  # num.bin :: Total number of bins to use in the scaled analysis (default is 20)
  
  # binning.method :: equal_width: each bin is of equal width along the x-axis [default]
  
  # mean.cutoff :: A two-length numeric vector with low- and high-cutoffs for feature means
  
  # dispersion.cutoff : A two-length numeric vector with low- and high-cutoffs for feature dispersions
  
  # nfeatures :: Number of genes to be selected as HVGs
  
  # Let's select seurat HVGs
  if (length(x = mean.cutoff) != 2 || length(x = dispersion.cutoff) != 2) {
    stop("Both 'mean.cutoff' and 'dispersion.cutoff' must be two numbers")
  }
  
  if (selection.method == "vst") {
    data <- GetAssayData(object = dataObject, slot = "counts")
    # if (ncol(x = data) < 1 || nrow(x = data) < 1) {
    if (IsMatrixEmpty(x = data)) {
      warning("selection.method set to 'vst' but count slot is empty; will use data slot instead")
      data <- GetAssayData(object = dataObject, slot = "data")
    }
  } else {
    data <- GetAssayData(object = dataObject, slot = "data")
  }
  
  hvf.info <- seurat_FindVariableFeatures.default(
    object = data,
    selection.method = selection.method,
    loess.span = loess.span,
    clip.max = clip.max,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    verbose = verbose
  )
  
  hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 0), ]
  
  if (selection.method == "vst") {
    hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
  } else {
    hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), , drop = FALSE]
  }
  
  selection.method <- switch(
    EXPR = selection.method,
    'mvp' = 'mean.var.plot',
    'disp' = 'dispersion',
    selection.method
  )
  
  top.features <- switch(
    EXPR = selection.method,
    'mean.var.plot' = {
      means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
      dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
      rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    },
    'dispersion' = head(x = rownames(x = hvf.info), n = nfeatures),
    'vst' = head(x = rownames(x = hvf.info), n = nfeatures),
    stop("Unknown selection method: ", selection.method)
  )
  
  return(top.features)
}




