as.GNARnet <- function(x){
  if(is.igraph(x)){
    return(igraphtoGNAR(x))
  }else{
    if(is.matrix(x)){
      return(matrixtoGNAR(x))
    }else{
      results <- vector(mode="logical", length=8)
      results[1] <- is.list(x)
      results[2] <- sum(names(x) == c("edges", "dist"))==2
      results[3] <- length(x$edges)==length(x$dist)
      if(all(results[1:3])){
        results[4] <- all(sapply(x$edges, length)==sapply(x$dist,length))
        tmp2 <- 1:length(x$edges)
        results[5] <- all(unlist(x$edges)%in%tmp2)
        results[6] <- all(sapply(x$edges, length)==sapply(lapply(x$edges, unique), length))
        results[7] <- all(unlist(x$dist)>0)
        selfloop <- function(n){
          return(n %in% x$edges[[n]])
        }
        results[8] <- all(sapply(tmp2, selfloop)==FALSE)
        if(all(results)){
          class(x) <- "GNARnet"
          return(x)
        }else{
          cat("Unable to convert this object")
        }
      }
    }
  }
}
