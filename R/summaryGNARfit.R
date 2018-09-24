summary.GNARfit <- function(object, ...){
  GNARobj <- object
  print(summary(GNARobj$mod))
  cat("GNAR BIC:", lmToBIC(GNARobj))
}
