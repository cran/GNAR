\name{matrixtoGNAR}
\alias{matrixtoGNAR}
\title{Converts an adjacency matrix to GNAR network}
\usage{
matrixtoGNAR(input.mat)
}
\description{
Converts an adjacency matrix to the \code{GNARnet} form for use as an input to GNAR functions.
}
\arguments{
\item{input.mat}{an adjacency matrix whose dimension equals the number of nodes in the resulting network.}
}
\details{
The values in the \code{$dist} list are the reciprocal of the values from the weighted adjacency matrix. Any self-loops (diagonal entries) and negatively weighted edges are removed.
}
\value{
\code{matrixtoGNAR} returns a \code{GNARnet} list with elements \code{$edges} and \code{$dist}.
}
\examples{
#Convert fiveNet to an adjacency matrix and back again
matrixtoGNAR(as.matrix(fiveNet))
}
