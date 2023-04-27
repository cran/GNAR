#) CHANGES START HERE
GNARXdesign <- function (vts = GNAR::fiveVTS, net = GNAR::fiveNet, alphaOrder = 2,
                          betaOrder = c(1,1), fact.var = NULL, globalalpha=TRUE,
                          tvnets=NULL, netsstart=NULL, lambdaOrder=NULL, xvts=NULL)
  #) CHANGES END HERE
  #) lambdaOrder is a length H vector containing the number of lags of the
  #) H exogenous regressors included in the GNARX equation.
  #) xvts is a length H list, where each element of the list is of same dimension
  #) as vts. If lambdaOrder remains null while xvts is not, xvts is ignored.
{
  #) CHANGES START HERE
  if(!is.null(lambdaOrder)){
    stopifnot(is.list(xvts))
    H <- length(lambdaOrder)
    stopifnot(H == length(xvts))
    stopifnot(floor(lambdaOrder) == lambdaOrder)
    stopifnot(min(lambdaOrder) >= 0)
    for(h in 1:H){
      stopifnot(all(dim(xvts[[h]]) == dim(vts)))
    }
  }
  #) CHANGES END HERE
  #) This ensures that each exogenous regressor is contained in a list,
  #) has a correponding lag order and
  #) all of the exogenous regressor data matrices have the same dimension as vts.
  #) Also ensure lag orders are nonnegative integers.
  stopifnot(is.GNARnet(net))
  stopifnot(ncol(vts) == length(net$edges))
  stopifnot(alphaOrder > 0)
  stopifnot(floor(alphaOrder) == alphaOrder)
  stopifnot(length(betaOrder) == alphaOrder)
  stopifnot(floor(betaOrder) == betaOrder)
  if(!is.null(fact.var)){
    stopifnot(length(fact.var) == length(net$edges))
    if(!globalalpha){
      stop("Use factors OR individual alphas")
    }
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
  #cat("Note: input net should contain distances (not weights)")
  #flip network so that NofNeighbours gives into node information
  netmat <- as.matrix(net, normalise=FALSE)
  if(!isSymmetric(netmat)){
    net <- as.GNARnet(t(netmat))
  }
  parNames <- parLoc <- NULL
  for (jj in 1:alphaOrder) {
    if(globalalpha){
      parNames <- c(parNames, paste("alpha", jj, sep = ""))
      parLoc <- c(parLoc, "a")
    }else{
      for(kk in 1:ncol(vts)){
        parNames <- c(parNames, paste("alpha", jj, "node", kk, sep=""))
        parLoc <- c(parLoc, "a")
      }
    }
    if (betaOrder[jj] > 0) {
      for (kk in 1:betaOrder[jj]) {
        parNames <- c(parNames, paste("beta", jj, ".",
                                      kk, sep = ""))
        parLoc <- c(parLoc, "b")
      }
    }
  }
  maxOrder <- alphaOrder
  #) CHANGES START HERE
  if(!is.null(lambdaOrder)){
    for(jj in 1:H){
      for(ii in 0:lambdaOrder[jj]){
        parNames <- c(parNames, paste("lambda", jj,".", ii, sep = ""))
        parLoc <- c(parLoc, paste("l", jj, sep = "")) 
      }
    }
    maxOrder <- max(maxOrder, max(lambdaOrder))
  }
  #) CHANGES END HERE
  #) Added a loop to assign names to the lambda coefficients of the exogenous
  #) regressors. Also changed the maxOrder just in case the lambda order
  #) exceeds the alpha order.
  
  predt <- nrow(vts) - maxOrder
  nnodes <- ncol(vts)
  
  #) CHANGES START HERE
  dmat <- matrix(0, nrow = predt * nnodes, ncol = length(parLoc), 
                  dimnames = list(NULL, parNames))
  #) CHANGES END HERE
  #) Changed the ncol argument to allow for exogenous regressor data.
  
  for (ii in 1:nnodes) {
    for (aa in 1:alphaOrder) {
      if(globalalpha){
        alphaLoc <- which(parLoc == "a")[aa]
      }else{
        alphaLoc <- which(parLoc=="a")[nnodes*(aa-1)+ii]
      }
      dmat[((predt * (ii - 1) + 1):(predt * ii)), alphaLoc] <- vts[((maxOrder +
                                                                       1 - aa):(predt + (maxOrder - aa))), ii]
    }
  }
  if (sum(betaOrder) > 0) {
    betaN <- NULL
    betaTimes <- rep(1:alphaOrder, betaOrder)
    for (jj in 1:alphaOrder) {
      #betaTimes <- c(betaTimes, rep(jj, betaOrder[jj]))
      if (betaOrder[jj] > 0) {
        betaN <- c(betaN, 1:betaOrder[jj])
      }
    }
    for (ii in 1:nnodes) {
      NofNei <- NofNeighbours(node = ii, stage = max(betaOrder),
                              net = net)
      Nei <- NofNei$edges
      Wei <- NofNei$dist
      if ((!is.null(Nei)) & (length(Nei) > 0)) {
        if (!is.null(Nei[[1]])&!is.na(Nei[[1]][1])) {
          Wei <- lapply(Wei, function(x){1/(x*sum(1/x))})
          # for (jj in 1:length(Wei)) {
          #   inv.dist <- 1/Wei[[jj]]
          #   Wei[[jj]] <- inv.dist/sum(inv.dist)
          # }
          for (bb in 1:sum(betaOrder)) {
            betaLoc <- which(parLoc == "b")[bb]
            if (length(Nei[[betaN[bb]]]) > 1) {
              # print(paste("node", ii, "betaN[bb]", betaN[bb]))
              # print("In length(Nei[[betaN[bb]]]) > 1")
              vts.cut <- vts[((maxOrder + 1 - betaTimes[bb]):(predt +
                                                                (maxOrder - betaTimes[bb]))), Nei[[betaN[bb]]]]
              for (kk in 1:nrow(vts.cut)) {
                if (any(is.na(vts.cut[kk, ]))) {
                  if (all(is.na(vts.cut[kk, ]))) {
                    #if there are no neighbours left at any time point, set to zero
                    vts.cut[kk, ] <- 0
                  }
                  else {
                    new.wei <- Wei[[betaN[bb]]][which(!is.na(vts.cut[kk,
                                                                     ]))]
                    new.wei <- new.wei/sum(new.wei)
                    sub.val <- vts.cut[kk, which(!is.na(vts.cut[kk,
                                                                ]))] %*% new.wei
                    vts.cut[kk, which(is.na(vts.cut[kk,
                                                    ]))] <- sub.val
                  }
                }
              }
              dmat[((predt * (ii - 1) + 1):(predt *
                                              ii)), betaLoc] <- vts.cut %*% Wei[[betaN[bb]]]
            }
            else {
              # print(paste("node", ii, "betaN[bb]", betaN[bb]))
              # print("In length(Nei[[betaN[bb]]]) > 1 else")
              if ((length(Nei[[betaN[bb]]]) == 1) &
                  (!is.na(Nei[[betaN[bb]]]))) {
                # print("In (length(Nei[[betaN[bb]]]) == 1) & (!is.na(Nei[[betaN[bb]]]))")
                vts.cut <- vts[((maxOrder +
                                   1 - betaTimes[bb]):(predt + (maxOrder -
                                                                  betaTimes[bb]))), Nei[[betaN[bb]]]]
                #and if this is missing at any time point, set to zero
                vts.cut[is.na(vts.cut)] <- 0
                dmat[((predt * (ii - 1) + 1):(predt *
                                                ii)), betaLoc] <- vts.cut * Wei[[betaN[bb]]]
              }
              else {
                dmat[((predt * (ii - 1) + 1):(predt *
                                                ii)), betaLoc] <- 0
              }
            }
          }
        }
        else {
          for (bb in 1:sum(betaOrder)) {
            betaLoc <- which(parLoc == "b")[bb]
            dmat[((predt * (ii - 1) + 1):(predt * ii)),
                 betaLoc] <- 0
          }
        }
      }
      else {
        for (bb in 1:sum(betaOrder)) {
          betaLoc <- which(parLoc == "b")[bb]
          dmat[((predt * (ii - 1) + 1):(predt * ii)),
               betaLoc] <- 0
        }
      }
    }
  }
  #) CHANGES START HERE
  if(!is.null(lambdaOrder)){
    assignCol <- ncol(dmat) - sum(lambdaOrder) - H + 1
    for(ii in 1:H){
      for(jj in 0:lambdaOrder[ii]){
        if(is.null(dim(xvts[[ii]]))){
          #) This if statement ensures that, if xvts[[ii]] is only a vector rather
          #) than a matrix, we don't use the vec() function, which throws up an 
          #) error when applied to a vector.
          dmat[, assignCol] <- matrix(xvts[[ii]][(maxOrder + 1 - jj):(nrow(xvts[[ii]]) - jj), ])
        }else{
          dmat[, assignCol] <- vec(xvts[[ii]][(maxOrder + 1 - jj):(nrow(xvts[[ii]]) - jj), ])
        }
        assignCol <- assignCol + 1
      }
    }
  }
  #) CHANGES END HERE
  #) This block of code incorporates the exogenous regressor data into the
  #) design matrix.
  
  
  if (is.null(fact.var)) {
    return(dmat)
  }else {
    #allow more than just two factors
    facun <- unique(fact.var)
    if(length(facun)==1){
      return(dmat)
    }else{
      dmcol <- ncol(dmat)
      dmatex <- dmat
      exnames <- paste(colnames(dmat), " '",facun[1],"'", sep="")
      for(ii in 2:length(facun)){ #duplicate matrix columns
        dmatex <- cbind(dmatex, dmat)
        #change names to reflect factors
        exnames <- c(exnames, paste(colnames(dmat), " '",facun[ii], "'", sep=""))
      }
      
      #for each unique factor, set other entries to 0
      
      for(ii in 1:length(facun)){
        dmatex[fact.var != facun[ii], ((ii-1)* dmcol + (1:dmcol))] <- 0
      }
      
      
      colnames(dmatex) <- exnames
      
      
      
      return(dmatex)
    }
    
  }
}
