residuals.GNARfit <- function(object,...){
  stopifnot(is.GNARfit(object))
  return(residToMat(object, nnodes=object$frbic$nnodes)$resid)
}
