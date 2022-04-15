## zscore function 
zscore <- function(input.data = input.data){
  input.data.rowsums <- rowSums(input.data)
  input.data.mean <- rowMeans(input.data)
  input.data.sd <- matrixStats::rowSds(as.matrix(input.data))
  names(input.data.sd) <- rownames(input.data)
  zscore <- (input.data-input.data.mean)/input.data.sd
  return(zscore)
}

zscore_result <- zscore(input.data  = input.data)
