\name{pnacf}
\alias{pnacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes the Partial Network Autocorrelation Function (PNACF)
%%  ~~function to do ... ~~
}
\description{Computes the PNACF for a choice of lag \eqn{h} and r-stage depth \eqn{r}, the PNACF is given by
\eqn{
  \mathrm{pnacf}(h, r) = \frac{\sum_{t=1}^{T - h} ( \boldsymbol{\hat{u}}_{t + h} - \boldsymbol{\overline{u}})^{T} \big ( \mathbf{W} \odot \mathbf{S}_r + \mathbf{I_d} \big )
    ( \boldsymbol{\hat{u}}_{t} - \boldsymbol{\overline{u}})}
    {\sum_{t=1}^{T}  ( \boldsymbol{\hat{u}}_{t} - \boldsymbol{\overline{u}})^{T} \big \{ \big (1 + \lambda \big) \mathbf{I_d} \big \} ( \boldsymbol{\hat{u}}_{t} - \boldsymbol{\overline{u}})},
}
where \eqn{\hat{\boldsymbol{X}}_{t}^{h - 1, r - 1} = \sum_{k = 1}^{h - 1} ( \hat{\alpha}_k \boldsymbol{X}_{t - k} + \sum_{s = 1}^{r - 1} \hat{\beta}_{ks} \boldsymbol{Z}_{t - k}^{s}  )}, \eqn{\boldsymbol{\hat{u}}_{t + h} = \boldsymbol{X}_{t + h} - \hat{\boldsymbol{X}}_{t + h}^{h - 1, r - 1}}, and \eqn{\boldsymbol{\hat{u}}_{t} = \boldsymbol{X}_{t} - \hat{\boldsymbol{X}}_{t}^{h - 1, r - 1}} are the empirical residuals corresponding to GNAR(h -1, [r-1, ..., r - 1]) fits, \eqn{\lambda} is the same as for the NACF; see \code{\link{nacf}}, and \eqn{\boldsymbol{\overline{u}}} is the mean of the fitted residuals
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  pnacf(h, s,  weight_matrix, stages_tensor, nts_data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{Lag (i.e., time-steps behind) at which the NACF is computed.}
  \item{s}{r-stage depth at which the NACF is computed (i.e., shortest distance between nodes).}
  \item{weight_matrix}{Weight matrix of the GNAR process, each entry corresponds to the weight between two nodes; see \code{\link{weights_matrix}} }
  \item{stages_tensor}{List of r-stage adjacency matrices \eqn{\mathbf{S}_r}, the order is ascending.}
  \item{nts_data}{Network time series observations, the number of rows is equal to the number of time steps, and the number of columns is equal to the number of series (variables).}
%%     ~~Describe \code{x} here~~
}
\value{If the network time series contains missing values, then the weights matrix and \eqn{\lambda} are adjusted, so that missing values do not contribute to the partial network autocorrelation. This is done by setting to zero the weights which correspond to a missing value and computing the new weight matrix and \eqn{\lambda} value.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Daniel Salnikov and Guy Nason
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{
Nason, G.P., Salnikov, D. and Cortina-Borja, M. (2023)
New tools for network time series with an application to COVID-19
hospitalisations. \url{https://arxiv.org/abs/2312.00530}
}

\examples{
#
# Compute the PNACF with respect to a stationary GNAR simulation
#
gnar_simulation <- GNARsim(n = 100, net=fiveNet, alphaParams = list(rep(0.35, 5)), 
                           betaParams = list(c(0.25)), sigma=1)
W = weights_matrix(fiveNet)
stages_list = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(fiveNet)), 3)
# PNACF
pnacf(3, 1, W, stages_list, gnar_simulation)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ pnacf }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
