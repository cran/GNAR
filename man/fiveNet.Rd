\name{fiveNet}
\alias{fiveNet}
\title{Example network}
\usage{
fiveNet
}
\description{
This example GNARnet object is a network containing five vertices connected by five undirected edges, each of length 1.
}
\format{
\code{fiveNet} is a GNARnet object containing \code{$edges} and \code{$dist}. \cr
\code{edges} is a list of length five, with \code{edges[[i]]} containing the vertices that node i is connected to. \cr
\code{dist} is a list of length five, with \code{dist[[i]]} containing the length of the vertices that node i is connected to.
}
\examples{
#Plot this network
data(fiveNode)
plot(fiveNet)
}
