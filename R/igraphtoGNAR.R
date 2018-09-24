igraphtoGNAR <- function(ig){
  stopifnot(is.igraph(ig))
  if(is.weighted(ig)){
    weimat <- get.adjacency(ig, attr="weight")
  }else{
    weimat <- get.adjacency(ig)
    weimat <- as.matrix(weimat)
  }
  net <- matrixtoGNAR(weimat)
  return(net)
}
