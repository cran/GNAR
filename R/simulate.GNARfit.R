simulate.GNARfit <- function(object, nsim=object$frbic$time.in, seed=NULL,
                             future=TRUE, set.noise=NULL, allcoefs=FALSE, ...){
  stopifnot(is.GNARfit(object))
  stopifnot(floor(nsim)==nsim)
  stopifnot(nsim>0)
  if(!is.null(seed)){
    stopifnot(floor(seed)==seed)
    set.seed(seed)
  }
  if(!is.null(object$frbic$fact.var)){
    stop("fact.var not currently supported with simulate")
  }
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  if(!is.null(set.noise)){
    sig <- set.noise
  }else{
    sig <- sigma(object$mod)
  }
  if(!allcoefs){
    nas <- is.na(object$mod$coefficients)
    pvs <- summary(object$mod)$coefficients[,4] < 0.05
    vals <- rep(0, length(pvs))
    vals[pvs] <- summary(object$mod)$coefficients[pvs,1]
    coefvec <- rep(0, length(nas))
    coefvec[(!nas)] <- vals

  }else{
    coefvec <- object$mod$coefficients
    coefvec[is.na(coefvec)] <- 0
  }

    #use GNARsim
    if(object$frbic$globalalpha){
      #global alpha has one alpha per time lag
      alphaout <-  vector(mode="list", length=object$frbic$alphas.in)
      betaout <- as.list(rep(0,length=object$frbic$alphas.in))
      count <- 1
      for(jj in 1:object$frbic$alphas.in){
        alphaout[[jj]] <- rep(coefvec[count], object$frbic$nnodes)
        if(object$frbic$betas.in[jj]>0){
          betaout[[jj]] <- coefvec[(count+1):(count+object$frbic$betas.in[jj])]
        }
        count <- count + object$frbic$betas.in[jj] + 1
      }

    }else{
      #multiple alphas per time lag
      alphaout <-  vector(mode="list", length=object$frbic$alphas.in)
      betaout <- as.list(rep(0,length=object$frbic$alphas.in))
      count <- 1
      for(jj in 1:object$frbic$alphas.in){
        alphaout[[jj]] <- coefvec[count:(count+object$frbic$nnodes-1)]
        if(object$frbic$betas.in[jj]>0){
          betaout[[jj]] <- coefvec[(count+object$frbic$nnodes):(count+
                                      object$frbic$nnodes+object$frbic$betas.in[jj]-1)]
        }
        count <- count + object$frbic$nnodes + object$frbic$betas.in[jj]
      }
    }
  if(!future){
    newser <- GNARsim(n=nsim, net = object$frbic$net.in,
                      alphaParams = alphaout, betaParams = betaout,
                      sigma=sig)
  }else{
    nnodes <- object$frbic$nnodes
    max.nei <- max(unlist(lapply(betaout, length)))
    nei.mats <- vector(mode="list", length=max.nei)
    #create weight matrices for neighbours
    for(ii in 1:max.nei){
      nei.mats[[ii]] <- as.matrix(x=object$frbic$net.in, stage=ii, normalise=TRUE)
      if(sum(nei.mats[[ii]])==0){
        warning("beta order too large for network, neighbour set ",ii," is empty")
      }
    }

    xx.init <- object$frbic$final.in
    ntimedep <- object$frbic$alphas.in
    stopifnot(nrow(xx.init)==ntimedep)

    xx.gen <- matrix(NA, nrow=nsim+ntimedep, ncol=nnodes)
    xx.gen[1:ntimedep,] <- xx.init

    for(tt in (ntimedep+1):(nsim+ntimedep)){

      for(ii in 1:ntimedep){
        if(ii==1){
          time.forecast <- alphaout[[ii]]*xx.gen[(tt-ii),]
        }else{
          tmpi <- alphaout[[ii]]*xx.gen[(tt-ii),]
          time.forecast <- time.forecast + tmpi
        }


      }

      nei.forecast <- 0
      beta.pos <- NULL
      for(aa in 1:ntimedep){
        bb <- length(betaout[[aa]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast <- nei.forecast + betaout[[aa]][dd]*xx.gen[tt-aa,]%*%t(nei.mats[[dd]])
          }
        }
      }
      xx.gen[tt,] <- time.forecast+nei.forecast+rnorm(n=object$frbic$nnodes, mean=0, sd=sig)
    }
    if(nsim==1){
      newser <- as.ts(t(xx.gen[(ntimedep+1):(nsim+ntimedep),]), start=1, end=nsim)
    }else{
      newser <- as.ts(xx.gen[(ntimedep+1):(nsim+ntimedep),], start=1, end=nsim)
    }
  }

  return(newser)
}
