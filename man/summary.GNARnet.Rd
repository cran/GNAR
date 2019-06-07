\name{summary.GNARnet}
\alias{summary.GNARnet}
\title{Summary function for GNAR networks}
\usage{
\method{summary}{GNARnet}(object, ...)
}
\description{
Prints brief information about a GNAR network.
}
\arguments{
\item{object}{the network\code{GNARnet} object associated with the time series, containing a list with entries \code{$edges} and \code{$dist}.}
\item{...}{additional arguments, unused here.}
}
\details{
S3 method for class "GNARnet".
}
\examples{
#print fiveNet summary information
summary(fiveNet)
}
