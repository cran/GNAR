\name{print.GNARnet}
\alias{print.GNARnet}
\title{Print function for GNAR networks}
\usage{
\method{print}{GNARnet}(x, ...)
}
\description{
Prints information about a GNAR network.
}
\arguments{
\item{x}{the network \code{GNARnet} object associated with the time series, containing a list with entries \code{$edges} and \code{$dist}.}
\item{...}{additional arguments, unused here.}
}
\details{
S3 method for class "GNARnet".
}
\examples{
#print fiveNet information
print(fiveNet)
}
