print.GNARnet <- function(x,...){
  stopifnot(is.GNARnet(x))
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  cat("GNARnet with", length(x$edges), "nodes", "\n")
  cat("edges:")
  edgect <- 0
  totedges <- length(unlist(x$edges))
  if(totedges<=100){
    for(ii in 1:length(x$edges)){
      if(!is.null(x$edges[[ii]])){
        for(jj in 1:length(x$edges[[ii]])){
          edgect <- edgect + 1
          cat(ii, "--", x$edges[[ii]][jj], " ", sep="")
          if(floor(edgect/10)==edgect/10){
            cat("\n")
            cat("     ")
          }
        }
      }

      #cat(paste(ii, x$edges[[ii]], sep="--", collapse=" ")," ")
    }
  }else{
    for(ii in 1:length(x$edges)){
      if(!is.null(x$edges[[ii]])){
        for(jj in 1:length(x$edges[[ii]])){
          if(edgect < 100){
            edgect <- edgect + 1
            cat(ii, "--", x$edges[[ii]][jj], " ", sep="")
            if(floor(edgect/10)==edgect/10){
              cat("\n")
              cat("     ")
            }
          }
        }
      }

      #cat(paste(ii, x$edges[[ii]], sep="--", collapse=" ")," ")
    }
    cat("... network too large to print all edges")

  }

  if(length(unique(unlist(x$dist)))==1){
    cat("\n", "edges of each of length ", unique(unlist(x$dist)), "\n")
  }else{
    cat("\n", "edges of unequal lengths", "\n")
  }
}
