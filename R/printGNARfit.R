print.GNARfit <- function(x, ...){
  stopifnot(is.GNARfit(x))
  cat("Model:", "\n")
  cat(paste("GNAR(", x$frbic$alphas.in,
            ",[", paste(x$frbic$betas.in,sep=",",
                        collapse=","),"])", sep=""),"\n")
  print(x$mod)
}
