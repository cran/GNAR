matrixtoGNAR <- function(input.mat){
  stopifnot(nrow(input.mat)==ncol(input.mat))
  input.mat <- as.matrix(input.mat)
  if(min(input.mat) < 0){
    cat("WARNING: negative entries present in original matrix, these will be removed\n")
    input.mat[input.mat<0] <- 0
  }

  if(sum(diag(input.mat))>0){
    cat("WARNING: diagonal entries present in original matrix, these will be removed\n")
    diag(input.mat) <- rep(0, nrow(input.mat))
  }
  edges <- dist <- vector(mode="list", length=nrow(input.mat))

  for(ii in 1:nrow(input.mat)){
    nz <- which(input.mat[ii,] != 0)
    if(length(nz)>0){
      edges[[ii]] <- nz
      dist[[ii]] <- 1/input.mat[ii,nz]
    }
  }
  out <- list(edges=edges, dist=dist)
  class(out) <- "GNARnet"
  return(out)
}
