\name{local_relevance_plot}
\alias{local_relevance_plot}
\alias{cross_correlation_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces a local neighbourhood relevance plot based on the distances in the underlying network.
%%  ~~function to do ... ~~
}
\description{Produces a local neighbourhood relevance plot based on the distances in the underlying network. The heat-map matrix should reflect clusters if a GNAR model is valid. The size of the clusters depends on the maximum r-stage depth for neighbourhood regression, as \eqn{r^*} gets larger, the clusters grow or intersect and cover more nodes. The relative strength of conditionally correlated nodes is \eqn{ \mathrm{rscc} (i, j) := \{ d(i,j) \}^{-1} \mathbb{I} \{ d(i, j) \leq r^* \} + \{2 d(i,j) \}^{-1} \mathbb{I} \{ r^* < d(i, j) \leq 2 r^* \}}. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  local_relevance_plot(network, r_star)
  cross_correlation_plot(h, vts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{GNAR network object, which is the underlying network for the time series under study.}
  \item{r_star}{Maximum active r-stage depth for neighbourhood regression.}
  \item{h}{The lag in the cross correlation plot.}
  \item{vts}{The vector time series to compute the cross correlation plot on.}
}
\value{Produces the local relevance plot. Does not return any values.
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
# Produces a local relevance plot, which is a heat-map matrix from a stationary 
# GNAR(1, [1]) simulation.
#
gnar_simulation <- GNARsim(n = 100, net=fiveNet, alphaParams = list(rep(0.35, 5)), 
        betaParams = list(c(0.25)), sigma=1)
# Active node plot
local_relevance_plot(fiveNet, 1)
# Compare to the cross-correlation plot at one-lag
cross_correlation_plot(1, gnar_simulation)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ sparse }
\keyword{ cross-correlation }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
