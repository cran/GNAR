is.GNARfit <- function(x){
  results <- vector(mode="logical", length=6)
  results[1] <- "mod" %in% names(x)
  results[2] <- any(c("y","ys") %in% names(x))
  results[3] <- any(c("dd","ds") %in% names(x))
  results[4] <- "frbic" %in% names(x)
  if(sum(results)==4){
    results[5] <- class(x$mod) == "lm"
    results[6] <- all(c("nnodes", "alphas.in", "betas.in",
                        "fact.var", "globalalpha") %in% names(x$frbic))
  }
  return(all(results))
}
