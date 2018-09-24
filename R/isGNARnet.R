is.GNARnet <- function(x){
  results <- vector(mode="logical", length=7)
  results[1] <- is.list(x)
  results[2] <- all(names(x) == c("edges", "dist"))
  results[3] <- length(x$edges)==length(x$dist)
  if(sum(results)==3){
    results[4] <- all(sapply(x$edges, length)==sapply(x$dist,length))
    tmp2 <- 1:length(x$edges)
    results[5] <- all(unlist(x$edges)%in%tmp2)
    results[6] <- class(x) == "GNARnet"
    results[7] <- all(sapply(x$edges, length)==sapply(lapply(x$edges, unique), length))
  }
  return(all(results))
}
