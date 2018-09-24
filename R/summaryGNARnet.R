summary.GNARnet <- function(object,...){
  stopifnot(is.GNARnet(object))
  cat("GNARnet with", length(object$edges), "nodes and", length(unlist(object$edges)),"edges")
  if(length(unique(unlist(object$dist)))==1){
    cat("\n", "of equal length ", unique(unlist(object$dist)))
  }else{
    cat("\n", "of unequal lengths")
  }
}
