\name{as.matrix.GNARnet}
\alias{as.matrix.GNARnet}
\title{Converts a GNAR networks into a weighted adjacency matrix}
\usage{
\method{as.matrix}{GNARnet}(x, stage=1, normalise=FALSE, ...)
}
\description{
Takes an input GNARnet and neighbour stage and outputs the corresponding adjacency matrix.
}
\arguments{
\item{x}{the network \code{GNARnet} object associated with the time series, containing a list with entries \code{$edges} and \code{$dist}.}
\item{stage}{the neighbour set that the adjacency matrix is created for.}
\item{normalise}{whether to normalise each to non-zero row to have sum one.}
\item{...}{additional arguments, unused here.}
}
\details{
S3 method for class "GNARnet".

With normalisation this is a non-invertible transform. See \link[GNAR]{NofNeighbours} for neighbour set definition.
}
\value{
\code{as.matrix} performed on a \code{GNARnet} returns a square matrix with the number of rows and columns equal to the length of the \code{$edges} list. Entry \code{i,j} of the matrix will be non-zero if node \code{j} is in the \code{stage} neighbour set of \code{i}.
}
\examples{
#fiveNet as an adjacency matrix
as.matrix(fiveNet)
}
