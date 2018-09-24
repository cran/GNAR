\name{plot.GNARnet}
\alias{plot.GNARnet}
\title{Plot function for GNAR networks}
\usage{
\method{plot}{GNARnet}(x, ...)
}
\description{
Plots a GNAR network using the 'igraph' package.
}
\arguments{
\item{x}{the network\code{GNARnet} object associated with the time series, containing a list with entries \code{$edges} and \code{$dist}.}
\item{...}{additional arguments for the \code{igraph} plotting function, see \link[igraph]{plot.igraph}.}
}
\details{
S3 method for class "GNARnet".
}
\examples{
#Plot fiveNet
data(fiveNode)
plot.GNARnet(fiveNet)
}
