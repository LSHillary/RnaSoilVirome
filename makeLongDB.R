makeLongDB <- function(path, pattern = NULL, chopFromFileName, sep = NULL, header){
  
  if(!require(data.table)){install.packages('data.table'); require(data.table)}
  
  # Get file names and locations
  files <- list.files(path = path, pattern = pattern)
  fileLocation <- paste(path, files, sep = "")
  
  # Make sure all files being read are larger than 0 bytes
  files <- files[file.size(fileLocation) != 0]
  fileLocation <- fileLocation[file.size(fileLocation) != 0]
  
  # Read files into a list. Standard seperator is tab, but others can be defined
  if(is.null(sep)){
    fileList <- lapply(fileLocation, FUN = function(x){
      data.frame(fread(file = x, sep = '\t', header = header))}
    )
  } else{
    fileList <- lapply(fileLocation, FUN = function(x){
      data.frame(fread(file = x, sep = sep, header = header))}
    )
  }
  
  # Add treatment names to the assemblies
  names(fileList) <- files
  
  # Put data from all assemblies into a data frame
  longDataFrame <- do.call("rbind", fileList)
  
  # Add a Sample column
  Sample <- sub(chopFromFileName,"",rep(names(fileList), sapply(fileList, nrow)))
  sampleDataFrame <- cbind(Sample, longDataFrame)
  return(sampleDataFrame)
}

RemoveContaminants <- function(X){
  IDX <- TableMapped[TableMapped$Sample == 'Negative Control',]
  IDX <- IDX[IDX$Covered_percent < X,]$vOTU
  dtX <- TableMapped[TableMapped$vOTU %in% IDX,]
  dtX <- dtX[dtX$Sample != 'Negative Control',]
  return(dtX)
}

IdentifyPresent <- function(X, J){
  J$Presence <- ifelse(J$Covered_percent < X, 0, 1)
  K <- data.table(J$vOTU, J$Sample, J$Presence)
  names(K) <- c("vOTU", "Sample", "Presence")
  K <- spread(K, key = Sample, value = Presence)
  L <- as.matrix(K[,-1], rownames = K$vOTU)
  L <- L[rowSums(L)>0,]
  return(L)
}