\name{fitted.GNARfit}
\alias{fitted.GNARfit}
\title{Function to return fitted values of GNARfit objects}
\usage{
\method{fitted}{GNARfit}(object,...)
}
\description{
\code{fitted.GNARfit} returns the fitted values of a GNARfit object as a matrix.
}
\arguments{
\item{object}{the output of a GNARfit or GNARpredict call}
\item{...}{additional arguments, unused here.}
}
\details{
S3 method for class "GNARfit".
}
\value{
\code{fitted.GNARfit} returns a \link[stats]{ts} object of fitted values, with \code{t-alphaOrder} rows and \code{nnodes} columns.
}
\examples{
#get the fitted values of the fiveNode GNAR fit
data(fiveNode)
fitted(GNARfit())
}
