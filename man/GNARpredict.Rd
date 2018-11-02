\name{GNARpredict}
\alias{GNARpredict}
\title{Fits and predicts using the GNAR model}
\usage{
GNARpredict(vts=GNAR::fiveVTS, net=GNAR::fiveNet, alphaOrder=2, betaOrder=c(1,1),
 fact.var=NULL, globalalpha=TRUE, tvnets=NULL, netsstart=NULL, ErrorIfNoNei=TRUE)
}
\description{
Fits the GNAR model to up to observation \code{t-1} and produces a within-sample prediction for time \code{t}.
}
\arguments{
\item{vts}{a matrix containing the multivariate time series to be modelled. The \code{i,j} entry of this matrix should be for time \code{i} and vertex/node \code{j}.}
\item{net}{the (first) network associated with the time series, containing a list with entries \code{$edges} and \code{$dist}. This network should have the same number of nodes as the number of columns of the \code{vts} matrix.}
\item{alphaOrder}{a non-negative integer specifying the maximum time-lag to model.}
\item{betaOrder}{a vector of length \code{alphaOrder} specifying the maximum neighbour set to model at each of the time-lags.}
\item{fact.var}{a vector of factors indicating which nodes belong to differents set with different parameters to be fitted.}
\item{globalalpha}{a TRUE/FALSE value indivating whether to use global alpha parameters.}
\item{tvnets}{a list of additional networks. Currently only NULL (the static network case) is supported.}
\item{netsstart}{a vector of times corresponding to the first time points for each network of \code{tvnets}. Currently only NULL (the static network case) is supported.}
\item{ErrorIfNoNei}{a TRUE/FALSE value indicating whether to stop the function call with an error if betaOrder specifies more neighbour sets than exist in the network. If FALSE the function will continue and some parameters will be NA.}
}
\details{
See \link[GNAR]{GNARfit} for GNAR model information. Note that the prediction is for the final time observation of the input data, to predict out-of-sample add a row of zeros to the input data matrix. Only coefficients with p-value smaller than 0.05 are used to calculate prediction. See \link[GNAR]{is.GNARnet} for \code{GNARnet} object information and example construction.
}
\value{
\item{pred}{the prediction for time \code{t}.}
\item{mod}{the \code{lm} output from fitting the GNAR model up to \code{t-1}.}
\item{ys}{the original response values up to \code{t-1}, with NAs left in.}
\item{ds}{the output of \code{GNARdesign} containing the design matrix up to \code{t-1}, with NAs left in.}
\item{ypred}{the original response values at \code{t}.}
\item{dpred}{the time \code{t} entries of the design matrix.}
\item{frbic}{inputs to the \code{GNARpredict} function.}
}
\examples{
#Fit and predict the fiveVTS data with the GNAR(2,[1,1]) model
data(fiveNode)
GNARpredict()
}
