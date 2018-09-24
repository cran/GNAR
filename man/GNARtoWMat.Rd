\name{GNARtoWMat}
\alias{GNARtoWMat}
\title{Converts a GNAR network to a weighted adjacency matrix}
\usage{
GNARtoWMat(net=GNAR::fiveNet, stage=1, normalise=FALSE)
}
\description{
Takes an input GNARnet and neighbour stage and outputs the corresponding adjacency matrix.
}
\arguments{
\item{net}{a GNARnet object containing \code{$edges} and \code{$dist}.}
\item{stage}{the neighbour set that the adjacency matrix is created for.}
\item{normalise}{whether to normalise each to non-zero row to have sum one.}
}
\details{
With normalisation this is a non-invertible transform. See \link[GNAR]{NofNeighbours} for neighbour set definition.
}
\value{
\code{GNARtoWMat} returns a square matrix with the number of rows and columns equal to the length of the \code{$edges} list. Entry \code{i,j} of the matrix will be non-zero if node \code{j} is in the \code{stage} neighbour set of \code{i}.
}
\examples{
#fiveNet as an adjacency matrix
data(fiveNode)
GNARtoWMat()
}
