GNARtoWMat <- function(net=GNAR::fiveNet, stage=1, normalise=FALSE){
  stopifnot(is.GNARnet(net))
  tmp.mat <- matrix(0, nrow=length(net$edges), ncol=length(net$edges))
  for(ii in 1:length(net$edges)){
    if(!is.null(net$edges[[ii]])){
      nnei <- NofNeighbours(node=ii, stage=stage, net=net)
      nwei <- 1/nnei$dist[[stage]]
      if(normalise){
        nwei <- nwei/sum(nwei)
      }
      tmp.mat[ii,nnei$edges[[stage]]] <- nwei
    }
  }
  return(tmp.mat)
}
