lmToBIC <- function(GNARobj.in=GNARpredict()){
  stopifnot(is.GNARfit(GNARobj.in))
  nnodes.in <- GNARobj.in$frbic$nnodes
  alphas.in <- GNARobj.in$frbic$alphas.in
  betas.in <- GNARobj.in$frbic$betas.in
  fact.var <- GNARobj.in$frbic$fact.var
  globalalpha <- GNARobj.in$frbic$globalalpha
  if(!is.null(fact.var)){
    cat("Factors not currently supported")
  }
  stopifnot(is.null(fact.var))
  # if(!globalalpha){
  #   cat("Individual alphas not currently supported")
  # }
  # stopifnot(globalalpha)
  stopifnot(is.logical(globalalpha))
  stopifnot(length(nnodes.in)==1)
  stopifnot(floor(nnodes.in)==nnodes.in)

  tmp.resid <- residToMat(GNARobj=GNARobj.in, nnodes=nnodes.in)$resid
  tmp.resid[is.na(tmp.resid)] <- 0 #shouldn't make a difference when comparing across p, s

  larg <- det((1/nrow(tmp.resid)) * t(tmp.resid) %*% tmp.resid)
  stopifnot(larg != 0)

  tmp1 <- log(larg)

  tot.time <- nrow(tmp.resid) + alphas.in
  stopifnot(tot.time != 0)

  if(globalalpha){
    tmp2 <- (  alphas.in + sum(betas.in) ) * log(tot.time) / tot.time
  }else{
    tmp2 <- ( ncol(tmp.resid) * alphas.in + sum(betas.in) ) * log(tot.time) / tot.time
  }

  return(tmp1 + tmp2)
}
