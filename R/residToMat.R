residToMat <- function(GNARobj=GNARfit(), nnodes=5){
  #allow either a fit or a predict output here
  #fit will have $y, predict will have $ys components
  stopifnot(floor(nnodes)==nnodes)
  stopifnot(is.GNARfit(GNARobj))
  if(is.null(GNARobj$ys)){
    yvec <- GNARobj$y
    dmat <- GNARobj$dd
  }else{
    yvec <- GNARobj$ys
    dmat <- GNARobj$ds
  }

  resid.gaps <- rep(0, length(yvec))
  resid.gaps[is.na(yvec)] <- 1
  resid.gaps[na.row(dmat)] <- 1

  resid.with.gaps <- rep(NA, length=length(yvec))
  resid.with.gaps[!resid.gaps] <- GNARobj$mod$residuals

  fit.with.gaps <- rep(NA, length=length(yvec))
  fit.with.gaps[!resid.gaps] <- GNARobj$mod$fitted

  resid.matrix <- matrix(resid.with.gaps,  ncol=nnodes, byrow=FALSE)
  fitted.matrix <- matrix(fit.with.gaps, ncol=nnodes, byrow=FALSE)

  return(list(resid=resid.matrix, fit=fitted.matrix))

}
