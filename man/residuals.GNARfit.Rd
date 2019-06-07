\name{residuals.GNARfit}
\alias{residuals.GNARfit}
\title{Function to return residuals of GNARfit objects}
\usage{
\method{residuals}{GNARfit}(object,...)
}
\description{
\code{residuals.GNARfit} returns the residuals of a GNARfit object as a matrix.
}
\arguments{
\item{object}{the output of a \link[GNAR]{GNARfit} call}
\item{...}{additional arguments, unused here.}
}
\details{
The function first checks if the object is of GNARfit class, then uses \link[GNAR]{residToMat} to return the residuals.
}
\value{
\code{residuals.GNARfit} returns a 'ts' object of residuals, with \code{t-alphaOrder} rows and \code{nnodes} columns.
}
\examples{
#get the residuals of the fiveNode GNAR fit
residuals(GNARfit())
}
