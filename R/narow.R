na.row <- function(mat){
  tmp <- apply(mat,1,function(x){any(is.na(x))})
  return(tmp)
}
