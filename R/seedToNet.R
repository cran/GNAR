seedToNet <- function(seed.no, nnodes=34, graph.prob=0.5){
  stopifnot(floor(nnodes)==nnodes)
  stopifnot(nnodes > 1)
  stopifnot(graph.prob <= 1)
  stopifnot(graph.prob >= 0)
  stopifnot(floor(seed.no)==seed.no)
  set.seed(seed.no)
  tmp.mat <- matrix(0, nrow=nnodes, ncol=nnodes)
  tmp.mat[lower.tri(tmp.mat)] <- rbinom(1,n=nnodes*(nnodes-1)/2, prob=graph.prob)
  tmp.mat[upper.tri(tmp.mat)] <- t(tmp.mat)[upper.tri(tmp.mat)]
  net.out <- matrixtoGNAR(tmp.mat)
  return(net.out)
}
