coef.GNARfit <- function(object,...){
  stopifnot(is.GNARfit(object))
  return(coef(object$mod))
}
