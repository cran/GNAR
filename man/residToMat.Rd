\name{residToMat}
\alias{residToMat}
\title{Converts the output of a GNARfit or GNARpredict to fitted and residual value matrices}
\usage{
residToMat(GNARobj=GNARpredict(), nnodes=5)
}
\description{
Unstacks the entries of the GNARfit or GNARpredict fitted and residual values to return matrices of a similar form to the multivariate time series input.
}
\arguments{
\item{GNARobj}{the output from the \link[GNAR]{GNARfit} or \link[GNAR]{GNARpredict} function}
\item{nnodes}{the number of nodes in the original network time series}
}
\details{
This function also replaces the NAs that were removed in fitting.
}
\value{
\item{resid}{is the matrix of residual values, with \code{t-alphaOrder} rows and \code{nnodes} columns.}
\item{fit}{is the matrix of fitted values, with \code{t-alphaOrder} rows and \code{nnodes} columns.}
}
\examples{
#Get residual and fitted matrices from GNARpredict fit of fiveVTS
data(fiveNode)
residToMat()
}
