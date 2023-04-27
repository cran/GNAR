#) CHANGES START HERE
GNARXfit <- function(vts=GNAR::fiveVTS, net=GNAR::fiveNet, alphaOrder=2, betaOrder=c(1,1), fact.var=NULL,
                    globalalpha=TRUE, tvnets=NULL, netsstart=NULL, ErrorIfNoNei=TRUE,
                    xvts=NULL, lambdaOrder=NULL){
  #) CHANGES END HERE
  #) Added arguments corresponding to exogenous regressor data (xvts) and lag orders
  #) (lambdOrder). xvts must be a list.
  #input checks
  
  #) CHANGES START HERE
  maxOrder <- alphaOrder
  if(!is.null(lambdaOrder)){
    stopifnot(is.list(xvts))
    H <- length(lambdaOrder)
    stopifnot(H == length(xvts))
    stopifnot(floor(lambdaOrder) == lambdaOrder)
    stopifnot(min(lambdaOrder) >= 0)
    for(h in 1:H){
      stopifnot(all(dim(xvts[[h]]) == dim(vts)))
    }
    maxOrder <- max(alphaOrder, lambdaOrder)
  }
  #) CHANGES END HERE
  #) This ensures that each exogenous regressor is contained in a list,
  #) has a correponding lag order and
  #) all of the exogenous regressor data matrices have the same dimension as vts.
  #) Also ensure lag orders are nonnegative integers.
  #) Also compute maxOrder, which replaces alphaOrder where appropriate in the below code.
  
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
  
  #) CHANGES START HERE
  frbic <- list(nnodes=length(net$edges),alphas.in=alphaOrder,betas.in=betaOrder,fact.var=fact.var,
                globalalpha=globalalpha, xtsp=tsp(vts), time.in=nrow(vts), net.in=net, final.in=vts[(nrow(vts)-maxOrder+1):nrow(vts),],
                lambdas.in = lambdaOrder)
  dmat <- GNARXdesign(vts=vts, net=net, alphaOrder=alphaOrder, betaOrder=betaOrder, fact.var=fact.var,
                     globalalpha=globalalpha, tvnets=tvnets, netsstart=netsstart,
                     xvts = xvts, lambdaOrder = lambdaOrder)
  #) CHANGES END HERE
  #) Replaced alphaOrder by maxOrder in the final.in element. Also included element lambdas.in.
  #) Extra two arguments included for GNARdesign. Here it's replaced by GNARXdesign,
  #) which can be renamed to GNARdesign once the new function has been approved.
  
  if(ErrorIfNoNei){
    if(any(apply(dmat==0, 2, all))){
      parname <- strsplit(names(which(apply(dmat==0, 2, all)))[1], split=NULL)[[1]]
      betastage <- parname[(which(parname==".")+1) :(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", betastage)
    }
  }
  
  #) CHANGES START HERE
  predt <- nrow(vts) - maxOrder
  #) CHANGES END HERE
  #) Replacing alphaOrder by maxOrder.
  
  yvec <- NULL
  for(ii in 1:length(net$edges)){
    #) CHANGES START HERE
    yvec <- c(yvec, vts[((maxOrder+1):(predt+maxOrder)),ii])
    #) CHANGES END HERE
    #) Replacing alphaOrder by maxOrder.
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
