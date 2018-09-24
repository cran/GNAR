\name{GNARdesign}
\alias{GNARdesign}
\title{Function to create the GNAR design matrix}
\usage{
GNARdesign(vts = GNAR::fiveVTS, net = GNAR::fiveNet, alphaOrder = 2, betaOrder = c(1,1),
 fact.var = NULL, globalalpha=TRUE, tvnets=NULL, netsstart=NULL)
}
\description{
Creates the design matrix necessary for fitting the GNAR model.
}
\arguments{
\item{vts}{a matrix or \link[stats]{ts} object containing the multivariate time series to be modelled. The \code{i,j} entry of this matrix should be for time \code{i} and vertex/node \code{j}.}
\item{net}{the (first) network associated with the time series, containing a list with entries \code{$edges} and \code{$dist}. This network should have the same number of nodes as the number of columns of the \code{vts} matrix.}
\item{alphaOrder}{a non-negative integer specifying the maximum time-lag to model.}
\item{betaOrder}{a vector of length \code{alphaOrder} specifying the maximum neighbour set to model at each of the time-lags.}
\item{fact.var}{a vector of factors indicating which nodes belong to each set with different parameters to be fitted.}
\item{globalalpha}{a TRUE/FALSE value indivating whether to use global alpha parameters.}
\item{tvnets}{a list of additional networks. Currently only NULL (the static network case) is supported.}
\item{netsstart}{a vector of times corresponding to the first time points for each network of \code{tvnets}. Currently only NULL (the static network case) is supported.}
}
\value{
\item{GNARdesign}{returns a matrix containing \code{(t-alphaOrder)nnodes} rows and a column for each parameter to be fitted. The columns are in time-lag order, eg for GNAR(2,[1,0]) the columns are \code{alpha1}, \code{beta1.1}, \code{alpha2}. When a factor variable is specified the columns are labelled with the factor.}
}
\examples{
#Design matrix to fit GNAR(2,[1,1]) to the fiveVTS data
data(fiveNode)
GNARdesign()
}
