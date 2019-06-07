coef.GNARfit <- function(object,...){
  stopifnot(is.GNARfit(object))
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  return(coef(object$mod))
}
