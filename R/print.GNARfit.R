print.GNARfit <- function(x, ...){
  stopifnot(is.GNARfit(x))
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  cat("Model:", "\n")
  cat(paste("GNAR(", x$frbic$alphas.in,
            ",[", paste(x$frbic$betas.in,sep=",",
                        collapse=","),"])", sep=""),"\n")
  print(x$mod)
}
