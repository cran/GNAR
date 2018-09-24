\name{is.GNARnet}
\alias{is.GNARnet}
\alias{as.GNARnet}
\title{Functions to check and create GNARnet objects}
\usage{
is.GNARnet(x)
as.GNARnet(x)
}
\description{
\code{is.GNARnet} returns either TRUE or FALSE according to a series of GNARnet checks. \code{as.GNARnet} returns a GNARnet object from an input weights matrix, 'igraph' object, or a GNARnet without assigned class.
}
\arguments{
\item{x}{the network to be tested or object to be converted}
}
\details{
The \code{is.GNARnet} function checks whether the network passes a series of tests:
\itemize{
  \item Is it a list containing $edges and $dist
  \item Are the $edges and $dist lists of the same length
  \item Are each of the elements of $edges the same length as the corresponding $dist element
  \item Do the edges only contain valid entries, 1,...,nnodes
  \item Is it labelled as \code{GNARnet} class
  \item Are no duplicate edges present
}
The \code{as.GNARnet} function converts \link[igraph]{igraph} objects to GNARnet form, other possible inputs are adjacency matrices, and lists with \code{$edges} and \code{$dist} entries of the correct form.
}
\value{
\code{is.GNARnet} returns \code{TRUE} or \code{FALSE} corresponding to passing the above tests.
\code{as.GNARnet} returns a \code{GNARnet} object.
}
\examples{
#check that the example network meets the criteria above
data(fiveNode)
is.GNARnet(fiveNet)

#convert to igraph and back again
as.GNARnet(GNARtoigraph(fiveNet))

#generate a new network with three nodes
#edges 1->2, 2->1, 2->3
#dist 1, 2, 1
#note 1->2 and 2->1 are of different lengths
threeEdge <- list(c(2), c(1,3), NULL)
threeDist <- list(c(1), c(2,1), NULL)
threeNet <- list(edges=threeEdge, dist=threeDist)
#check if this is a GNARnet
is.GNARnet(threeNet)
#use as.GNARnet to change the class
threeNet <- as.GNARnet(threeNet)
#check if this is a GNARnet now
is.GNARnet(threeNet)
}
