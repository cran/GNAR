GNARsim <- function(n=200, net=GNAR::fiveNet, alphaParams=list(c(rep(0.2,5))), betaParams=list(c(0.5)), sigma=1,
                    tvnets=NULL, netsstart=NULL){
  #use 0s in betaParams to discount dependencies
  stopifnot(is.GNARnet(net))
  stopifnot(length(alphaParams)==length(betaParams))
  stopifnot(!is.null(betaParams))
  stopifnot(length(alphaParams[[1]])==length(net$edges))

  if(!is.null(tvnets)){
    cat("Time-varying nets not currently supported")
  }
  stopifnot(is.null(tvnets))

  nnodes <- length(net$edges)
  max.nei <- max(unlist(lapply(betaParams, length)))
  nei.mats <- vector(mode="list", length=max.nei)
  #create weight matrices for neighbours
  for(ii in 1:max.nei){
    nei.mats[[ii]] <- as.matrix(x=net, stage=ii, normalise=TRUE)
    if(sum(nei.mats[[ii]])==0){
      warning("beta order too large for network, neighbour set ",ii," is empty")
    }
  }
  #print(length(alphaParams))

  #seed the process from normal dist with mean 0 and sigma as given
  #do this for as many time points as needed for alpha
  ntimedep <- length(alphaParams)
  #print(length(alphaParams))

  xx.init <- matrix(rnorm(nnodes*ntimedep, mean=0, sd=sigma), nrow=ntimedep, ncol=nnodes)

  xx.gen <- matrix(NA, nrow=n+50, ncol=nnodes)
  xx.gen[1:ntimedep,] <- xx.init
  # print(length(alphaParams))
  # print("just before the tt indexed loop")
  for(tt in (ntimedep+1):(n+50)){
    # print(tt)
    # print(alphaParams)
    # print(rev(alphaParams))
    for(ii in 1:ntimedep){
      if(ii==1){
        time.forecast <- alphaParams[[ii]]*xx.gen[(tt-ii),]
      }else{
        tmpi <- alphaParams[[ii]]*xx.gen[(tt-ii),]
        time.forecast <- time.forecast + tmpi
      }


    }

    nei.forecast <- 0
    beta.pos <- NULL
    for(aa in 1:ntimedep){
      bb <- length(betaParams[[aa]])
      if(bb>0){
        for(dd in 1:bb){
          nei.forecast <- nei.forecast + betaParams[[aa]][dd]*xx.gen[tt-aa,]%*%t(nei.mats[[dd]])
        }
      }
    }
    xx.gen[tt,] <- time.forecast+nei.forecast+rnorm(nnodes, mean=0, sd=sigma)
  }
  return(as.ts(xx.gen[51:(n+50),]))
}
