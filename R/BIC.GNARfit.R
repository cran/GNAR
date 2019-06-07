BIC.GNARfit <- function(object,...){
    stopifnot(is.GNARfit(object))
    nnodes.in <- object$frbic$nnodes
    alphas.in <- object$frbic$alphas.in
    betas.in <- object$frbic$betas.in
    fact.var <- object$frbic$fact.var
    tot.time <- object$frbic$time.in
    globalalpha <- object$frbic$globalalpha
    dotarg <- list(...)
    if(length(dotarg)!=0){
      if(!is.null(names(dotarg))){
        warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
      }else{
        warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
      }
    }
    if(!is.null(fact.var)){
      f.in <- length(unique(fact.var))
    }else{
      f.in <- 1
    }
    stopifnot(is.logical(globalalpha))
    stopifnot(length(nnodes.in)==1)
    stopifnot(floor(nnodes.in)==nnodes.in)
    stopifnot(tot.time != 0)

    tmp.resid <- residToMat(GNARobj=object, nnodes=nnodes.in)$resid
    tmp.resid[is.na(tmp.resid)] <- 0 #shouldn't make a difference when comparing across p, s

    larg <- det((1/tot.time) * t(tmp.resid) %*% tmp.resid)
    stopifnot(larg != 0)

    tmp1 <- log(larg)

    if(globalalpha){
      tmp2 <- f.in * (  alphas.in + sum(betas.in) ) * log(tot.time) / tot.time
    }else{
      tmp2 <- ( ncol(tmp.resid) * alphas.in + sum(betas.in) ) * log(tot.time) / tot.time
    }

    return(tmp1 + tmp2)
}
