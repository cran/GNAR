GNARpredict <- function(vts=GNAR::fiveVTS, net=GNAR::fiveNet, alphaOrder=2, betaOrder=c(1,1), fact.var=NULL,
                        globalalpha=TRUE, tvnets=NULL, netsstart=NULL, ErrorIfNoNei=TRUE){
  #the last row of vts will be predicted
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
  if(ErrorIfNoNei){
    if(any(apply(dmat==0, 2, all))){
      parname <- strsplit(names(which(apply(dmat==0, 2, all)))[1], split=NULL)[[1]]
      betastage <- parname[(which(parname==".")+1) :(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", betastage)
    }
  }

  predt <- nrow(vts)-alphaOrder
  yvec <- NULL
  for(ii in 1:length(net$edges)){
    yvec <- c(yvec, vts[((alphaOrder+1):(predt+alphaOrder)),ii])
  }
  # if(any(is.na(yvec))|any(is.na(dmat))){
  #   cat("\n")
  #   cat("Note: NAs present - lm removes these from model")
  # }

  #strip out final values from dmat and yvec
  dmat.pred <- dmat[predt*(1:length(net$edges)),]
  yvec.pred <- yvec[predt*(1:length(net$edges))]

  dmat.st <- dmat[-c(predt*(1:length(net$edges))),] #use .st to fit model, then .pred to get new value
  yvec.st <- yvec[-c(predt*(1:length(net$edges)))]

  if(is.vector(dmat.st)){
    dmat.st <- as.matrix(dmat.st, ncol=1)
  }
  if(is.vector(dmat.pred)){
    dmat.pred <- as.matrix(dmat.pred, ncol=1)
  }

  if(sum(is.na(yvec.st))>0){
    yvec2.st <- yvec.st[!is.na(yvec.st)]
    dmat2.st <- dmat.st[!is.na(yvec.st),]
    modNoIntercept <- lm(yvec2.st~dmat2.st+0)

  }else{
    modNoIntercept <- lm(yvec.st~dmat.st+0)
  }


  #use only significant parameters in fit
  coef.locs <- function(mod){
    #accomodate NAs in coefficients
    nas <- is.na(mod$coefficients)
    pvs <- summary(mod)$coefficients[,4] < 0.05
    if(all(!nas)){
      return(pvs)
    }else{
      out <- !nas
      out[out] <- pvs
      return(out)
    }
    
  }

  use.coef <- coef.locs(modNoIntercept)
  #dmat.pred <- dmat.pred

  if(is.vector(dmat.pred)){
    dmat.pred <- matrix(dmat.pred, ncol=1)
    pred <- t(dmat.pred[,use.coef]) * coef(modNoIntercept)[use.coef]

  }else{
    if(sum(use.coef)==1){
      pred <- t(dmat.pred[,use.coef]) * coef(modNoIntercept)[use.coef]

    }else{

      pred <- dmat.pred[,use.coef] %*% coef(modNoIntercept)[use.coef]
    }
  }
  out <- list(pred=c(pred), mod=modNoIntercept, ys=yvec.st, ds=dmat.st, ypred=yvec.pred, dpred=dmat.pred, frbic=frbic)
  class(out) <- "GNARfit"
  return(out)

}
