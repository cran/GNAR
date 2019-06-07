summary.GNARfit <- function(object, ...){
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  print(summary(object$mod))
  cat("GNAR BIC:", BIC.GNARfit(object))
}
