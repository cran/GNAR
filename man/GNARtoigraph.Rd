\name{GNARtoigraph}
\alias{GNARtoigraph}
\title{Converts a GNAR network to a weighted igraph object}
\usage{
GNARtoigraph(net=GNAR::fiveNet, stage=1, normalise=FALSE)
}
\description{
Takes an input network and neighbour stage and returns it in \link[igraph]{igraph} form.
}
\arguments{
\item{net}{a GNARnet object containing \code{$edges} and \code{dist}.}
\item{stage}{the neighbour set that the adjacency matrix is created for.}
\item{normalise}{whether to normalise each to non-zero row to have sum one.}
}
\details{
With normalisation this is a non-invertible transform. See \link[GNAR]{NofNeighbours} for neighbour set definition. See \link[GNAR]{is.GNARnet} for \code{GNARnet} object information and example construction.
}
\value{
\code{GNARtoigraph} returns an 'igraph' object with weights as the inverse distances of the input network.
}
\examples{
#fiveNet as an igraph object
data(fiveNode)
GNARtoigraph()
}
