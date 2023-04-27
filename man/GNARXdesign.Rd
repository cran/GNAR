\name{GNARXdesign}
\alias{GNARXdesign}
\title{Function to create the GNARX design matrix}
\usage{
GNARXdesign(vts = GNAR::fiveVTS, net = GNAR::fiveNet, alphaOrder = 2, betaOrder = c(1,1),
 fact.var = NULL, globalalpha=TRUE, tvnets=NULL, netsstart=NULL, lambdaOrder=NULL, 
 xvts=NULL)
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
\item{lambdaOrder}{a vector of the same length as \code{xvts} containing non-negative integers specifying the maximum time-lag of exogenous regressors to model. The \code{h}th element of the vector refers to the maximum lag of the \code{h}th exogenous regressor.}
\item{xvts}{a list of matrices containing values of the exogenous regressors for each vertex/node. The \code{i,j} entry of the \code{h}th element of the list refers to the value of the \code{h}th exogenous regressor for time \code{i} and vertex/node \code{j}.}
}
\value{
\item{GNARdesign}{returns a matrix containing \code{(t-alphaOrder)nnodes} rows and a column for each parameter to be fitted. The columns are in time-lag order, eg for GNAR(2,[1,0]) the columns are \code{alpha1}, \code{beta1.1}, \code{alpha2}. When a factor variable is specified the columns are labelled with the factor. If exogenous regressors are included, the matrix takes the form given in Equation 10 of the supplement to the referenced paper.}
}
\references{Nason G.P. and Wei J. (2022) Quantifying the economic response to COVID-19 mitigations and death rates via forecasting Purchasing Managers’ Indices using Generalised Network Autoregressive models with exogenous variables. \emph{Journal of the Royal Statistical Society Series A}, \bold{185}, 1778--1792.
}
\examples{

set.seed(1)
n = 1000
xvts=list()
xvts[[1]] = matrix(rnorm(5*n, mean=0, sd=2), nrow=n, ncol=5)
xvts[[2]] = matrix(rnorm(5*n, mean=0, sd=2), nrow=n, ncol=5)
lambdaParams=list()
lambdaParams[[1]] = c(0.5, -0.5)
lambdaParams[[2]] = c(0.3, 0.1)

# Simulate the GNARX using the exogenous variables xvts with associated parameters lambdaParams

Y_data <- GNARXsim(n=n, net=GNAR::fiveNet, alphaParams=list(c(rep(0.2,5))), betaParams=list(c(0.5)),
                      sigma=1, xvts=xvts, lambdaParams=lambdaParams)

#Design matrix to fit GNARX(2,[1,1]) to the fiveVTS data

Xdesign <- GNARXdesign(vts = Y_data, net = GNAR::fiveNet, globalalpha = TRUE, alphaOrder = 1,
                betaOrder = 1, xvts = xvts, lambdaOrder = c(1,1))

Xdesign
}
