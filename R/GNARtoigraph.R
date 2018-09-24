GNARtoigraph<- function(net=GNAR::fiveNet, stage=1, normalise=FALSE){
  stopifnot(is.GNARnet(net))
  weimat <- GNARtoWMat(net=net, stage=stage, normalise=normalise)
  #set as undirected igraph if adjacency matrix is symmetric
  if(all(isSymmetric(weimat))){
    tmp <- graph.adjacency(weimat, mode="undirected", weighted=TRUE)
  }else{
    tmp <- graph.adjacency(weimat, mode="directed", weighted=TRUE)
  }
  return(tmp)
}
