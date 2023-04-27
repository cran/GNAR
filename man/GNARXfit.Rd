\name{GNARXfit}
\alias{GNARXfit}
\title{Fitting function for GNARX models}
\usage{
GNARXfit(vts=GNAR::fiveVTS, net=GNAR::fiveNet, alphaOrder=2, betaOrder=c(1,1),
 fact.var=NULL, globalalpha=TRUE, tvnets=NULL, netsstart=NULL, ErrorIfNoNei=TRUE, 
 xvts=NULL, lambdaOrder=NULL)
}
\description{
Fits the GNARX model to the given inputs using \code{GNARdesignX} and \code{lm}.
}
\arguments{
\item{vts}{a matrix containing the multivariate time series to be modelled. The \code{i,j} entry of this matrix should be for time \code{i} and vertex/node \code{j}.}
\item{net}{the (first) network associated with the time series, containing a list with entries \code{$edges} and \code{$dist}. This network should have the same number of nodes as the number of columns of the \code{vts} matrix.}
\item{alphaOrder}{a non-negative integer specifying the maximum time-lag to model.}
\item{betaOrder}{a vector of length \code{alphaOrder} specifying the maximum neighbour set to model at each of the time-lags.}
\item{fact.var}{a vector of factors indicating which nodes belong to different sets with different parameters to be fitted.}
\item{globalalpha}{a TRUE/FALSE value indivating whether to use global alpha parameters.}
\item{tvnets}{a list of additional networks. Currently only NULL (the static network case) is supported.}
\item{netsstart}{a vector of times corresponding to the first time points for each network of \code{tvnets}. Currently only NULL (the static network case) is supported.}
\item{ErrorIfNoNei}{a TRUE/FALSE value indicating whether to stop the function call with an error if betaOrder specifies more neighbour sets than exist in the network. If FALSE the function will continue and some parameters will be NA.}
\item{xvts}{a list of matrices containing values of the exogenous regressors for each vertex/node. The \code{i,j} entry of the \code{h}th element of the list refers to the value of the \code{h}th exogenous regressor for time \code{i} and vertex/node \code{j}.}
\item{lambdaOrder}{a vector of the same length as \code{xvts} containing non-negative integers specifying the maximum time-lag of exogenous regressors to model. The \code{h}th element of the vector refers to the maximum lag of the \code{h}th exogenous regressor.}
}


\details{
The GNAR model of order \eqn{(p,S)}{(p,S)} is defined as
\deqn{X_{i,t} = \sum_{j=1}^p \left( \alpha_{i,j} X_{i,t-j} +\sum_{c=1}^C \sum_{r=1}^{S_j} \beta_{j,r, c}\sum_{q \in {N}^{(r)}_t(i)} \omega_{i,q,c}^{(t)} X_{q,t-j} \right) + u_{i,t}}{X_{i,t} = sum_{j=1}^p ( \alpha_{i,j} X_{i,t-j} +sum_{c=1}^C sum_{r=1}^{S_j} \beta_{j,r, c} sum_{q \in {N}^{(r)}_t(i)} \omega_{i,q,c}^{(t)} X_{q,t-j} ) + u_{i,t},}
where \eqn{p} is the maximum time lag, \eqn{S=(S_1,...,S_p)}{S=(S_1,...,S_p)} and \eqn{S_j } is the maximum stage of neighbour dependence for time lag \eqn{j}, \eqn{{N}^{(r)}_t(i)}{N^(r)_t (i)} is the \eqn{r}th stage neighbour set of node \eqn{i} at time \eqn{t}, \eqn{\omega_{i,q,c}^{(t)}}{\omega_{i,q,c}^(t)} is the connection weight between node \eqn{i} and node \eqn{q} at time \eqn{t} if the path corresponds
to covariate \eqn{c}. Here, we consider a sum from one to zero to be zero and \eqn{\{u_{i,t}\}}{{u_{i,t}}}, are assumed to be independent and identically distributed at each node \eqn{i}, with mean zero and variance \eqn{\sigma_i^2}{sigma_i^2}.
Currently, only a single network GNAR model can be fitted.
The connection weight, \eqn{\omega_{i,q,c}^{(t)}}{\omega_{i,q,c}^(t)}, is the inverse of the distance between nodes \code{i} and \code{q}, normalised so that they sum to 1 for each \code{i,t}.
See \link[GNAR]{is.GNARnet} for \code{GNARnet} object information and example construction.

A GNARX process of order \eqn{p}, neighbourhood order vector \eqn{s} of length \eqn{p} and \eqn{H} node-specific time series exogenous variables with order vector \eqn{p'}, denoted 
 GNARX \eqn{(p, s, p')}, is given by
 
\deqn{
    Y_{i,t}=\sum_{j=1}^{p}\left(\alpha_{i,j}Y_{i,t-j} + \sum_{r=1}^{s_j}  \beta_{j,r} \sum_{q\in\mathcal{N}^{(r)}(i)}\omega_{i,q}Y_{q,t-j}\right)+ {\color{blue} \sum_{h=1}^H
    \sum_{j'=0}^{p'_h}
    \lambda_{h,j'} X_{h,i,t-j'}+u_{i,t}},
}

 where \eqn{ u_{i, t}} are assumed to be set of mutually uncorrelated random variables with
 mean zero and variance of \eqn{\sigma^2}.



}
\value{
\item{mod}{the \code{lm} output from fitting the GNAR model.}
\item{y}{the original response values, with NAs left in.}
\item{dd}{the output of \code{GNARdesign} containing the design matrix, with NAs left in.}
\item{frbic}{inputs to other \code{GNAR} functions.}
}
\references{Knight, M.I., Nunes, M.A. and Nason, G.P. Modelling, detrending and decorrelation of network time series.
	\href{https://arxiv.org/abs/1603.03221}{arXiv preprint}.\cr\cr
Knight, M.I., Leeming, K., Nason, G.P. and Nunes, M. A. (2020) Generalised Network Autoregressive Processes and the GNAR package. \emph{Journal of Statistical Software}, \bold{96} (5), 1--36.\cr\cr
Nason G.P. and Wei J. (2022) Quantifying the economic response to COVID-19 mitigations and death rates via forecasting Purchasing Managers’ Indices using Generalised Network Autoregressive models with 
exogenous variables. \emph{Journal of the Royal Statistical Society Series A}, \bold{185}, 1778--1792.
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

# fit a GNARX to the model

model <- GNARXfit(vts = Y_data, net = GNAR::fiveNet,globalalpha = TRUE, alphaOrder = 1,
                betaOrder = 1, xvts = xvts, lambdaOrder = c(1,1))

# look at the residuals
residToMat(model)$resid
}
