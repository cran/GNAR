summary.GNARnet <- function(object,...){
  stopifnot(is.GNARnet(object))
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  cat("GNARnet with", length(object$edges), "nodes and", length(unlist(object$edges)),"edges")
  if(length(unique(unlist(object$dist)))==1){
    cat("\n", "of equal length ", unique(unlist(object$dist)))
  }else{
    cat("\n", "of unequal lengths")
  }
}
