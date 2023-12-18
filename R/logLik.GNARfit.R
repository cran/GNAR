logLik.GNARfit <- function(object,...){
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  t <- object$frbic$time.in
  n <- object$frbic$nnodes
  eps <- residuals(object)
  eps[is.na(eps)] <- 0
  Sigma <- crossprod(eps)/t
  stopifnot(floor(t)==t)
  stopifnot(floor(n)==n)
  ll <- -(t*n) * log(2*pi)/2 - (t/2)*log(det(Sigma)) - sum(diag(eps %*% solve(Sigma) %*% t(eps)))/2
  class(ll) <- "logLik"
  attr(ll, "df") <- length(coef(object))
  return(ll)
}


