predict.GNARfit <- function(object, n.ahead=1, ...){
  stopifnot(is.GNARfit(object))
  if(!is.null(object$frbic$fact.var)){
    stop("fact.var not currently supported with predict")
  }
  return(simulate.GNARfit(object, nsim=n.ahead, future=TRUE, set.noise=0, ...))
}
