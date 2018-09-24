GNARfit <- function(vts=GNAR::fiveVTS, net=GNAR::fiveNet, alphaOrder=2, betaOrder=c(1,1), fact.var=NULL,
                    globalalpha=TRUE, tvnets=NULL, netsstart=NULL){
  #input checks
  stopifnot(is.GNARnet(net))
  stopifnot(ncol(vts) == length(net$edges))
  stopifnot(alphaOrder > 0)
  stopifnot(floor(alphaOrder) == alphaOrder)
  stopifnot(length(betaOrder) == alphaOrder)
  stopifnot(floor(betaOrder) == betaOrder)
  if(!is.null(fact.var)){
    stopifnot(length(fact.var) == length(net$edges))
    # if(sum(fact.var %in% c(0,1))!=length(fact.var)){
    #   cat("More than two (0/1) factor variables not yet supported")
    # }
    # stopifnot(sum(fact.var %in% c(0,1))==length(fact.var))
  }
  stopifnot(is.matrix(vts))
  # if(!globalalpha){
  #   cat("Individual alphas not yet supported")
  # }
  # stopifnot(globalalpha)
  stopifnot(is.logical(globalalpha))
  if(!is.null(tvnets)){
    cat("Time-varying networks not yet supported")
  }
  stopifnot(is.null(tvnets))
  useNofNei <- 1
  #cat("Note: input net should contain distances (not weights)")
  #end of input checks
  frbic <- list(nnodes=length(net$edges),alphas.in=alphaOrder,betas.in=betaOrder,fact.var=fact.var,
                globalalpha=globalalpha, xtsp=tsp(vts))
  dmat <- GNARdesign(vts=vts, net=net, alphaOrder=alphaOrder, betaOrder=betaOrder, fact.var=fact.var,
                     globalalpha=globalalpha, tvnets=tvnets, netsstart=netsstart)
  predt <- nrow(vts)-alphaOrder
  yvec <- NULL
  for(ii in 1:length(net$edges)){
    yvec <- c(yvec, vts[((alphaOrder+1):(predt+alphaOrder)),ii])
  }
  # if(any(is.na(yvec))|any(is.na(dmat))){
  #   cat("\\n")
  #   cat("Note: NAs present - lm removes these from model")
  # }
  if(sum(is.na(yvec))>0){
    yvec2 <- yvec[!is.na(yvec)]
    dmat2 <- dmat[!is.na(yvec),]
    modNoIntercept <- lm(yvec2~dmat2+0)

  }else{
    modNoIntercept <- lm(yvec~dmat+0)
  }
  out <- list(mod=modNoIntercept, y=yvec, dd=dmat, frbic=frbic)
  class(out) <- "GNARfit"
  return(out)
}
