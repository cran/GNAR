plot.GNARnet <- function(x, ...){
  ignet <- GNARtoigraph(net=x)
  plot(ignet,...)
}
