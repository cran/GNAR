\name{GNARXsim}
\alias{GNARXsim}
\title{Simulates a GNARX process}
\usage{
GNARXsim(n=200, net=GNAR::fiveNet, alphaParams=list(c(rep(0.2,5))),
	betaParams=list(c(0.5)), sigma=1, tvnets=NULL, netsstart=NULL, xvts=NULL, 
	lambdaParams=NULL)
}
\description{
Simulates a GNAR process with Normally distributed innovations.
}
\arguments{
\item{n}{time length of simulation.}
\item{net}{network used for the GNAR simulation.}
\item{alphaParams}{a list containing vectors of auto-regression parameters for each time-lag.}
\item{betaParams}{a list of equal length as \code{alphaParams} containing the network-regression parameters for each time-lag.}
\item{sigma}{the standard deviation for the innovations.}
\item{tvnets}{a list of additional networks. Currently only NULL (the static network case) is supported.}
\item{netsstart}{a vector of times corresponding to the first time points for each network of \code{tvnets}. 
Currently only NULL (the static network case) is supported.}
\item{xvts}{a list of matrices containing values of the exogenous regressors for each vertex/node. 
The \code{i,j} entry of the \code{h}th element of the list refers to the value of the \code{h}th exogenous 
regressor for time \code{i} and vertex/node \code{j}.}
\item{lambdaParams}{a list containing vectors of parameters associated to effect of the exogenous regressor variables for each time-lag.}
}
\details{
Parameter lists should not be NULL, set unused parameters to be zero. See \link[GNAR]{GNARXfit} for model description.
}
\value{
\code{GNARXsim} returns the multivariate time series as a \link[stats]{ts} object, with \code{n} rows and a column for each of the nodes in the network.
}
\references{Knight, M.I., Nunes, M.A. and Nason, G.P. Modelling, detrending and decorrelation of network time series.
	\href{https://arxiv.org/abs/1603.03221}{arXiv preprint}.\cr\cr
Knight, M.I., Leeming, K., Nason, G.P. and Nunes, M. A. (2020) Generalised Network Autoregressive Processes and the GNAR package. \emph{Journal of Statistical Software}, \bold{96} (5), 1--36.\cr\cr
Nason G.P. and Wei J. (2022) Quantifying the economic response to COVID-19 mitigations and death rates via forecasting Purchasing Managers’ Indices using Generalised Network Autoregressive models with
exogenous variables. \emph{Journal of the Royal Statistical Society Series A}, \bold{185}, 1778--1792.
}
\examples{

#Simulate a GNARX process with the fiveNet network

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

# now try to refit the model
 
model <- GNARXfit(vts = Y_data, net = GNAR::fiveNet,globalalpha = TRUE, alphaOrder = 1, 
		betaOrder = 1, xvts = xvts, lambdaOrder = c(1,1))

model
}
