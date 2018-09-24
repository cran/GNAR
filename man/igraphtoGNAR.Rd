\name{igraphtoGNAR}
\alias{igraphtoGNAR}
\title{Converts an igraph to GNAR network}
\usage{
igraphtoGNAR(ig)
}
\description{
Converts an 'igraph' to the \code{GNARnet} form for use as an input to GNAR functions.
}
\arguments{
\item{ig}{an 'igraph' object.}
}
\details{
The values in the \code{$dist} list are the reciprocal of the values from the weighted adjacency matrix.
}
\value{
\code{igraphtoGNAR} returns a \code{GNARnet}: a list with elements \code{$edges} and \code{$dist}.
}
\examples{
#Convert fiveNet to igraph and back again
data(fiveNode)
igraphtoGNAR(GNARtoigraph(fiveNet))
}
