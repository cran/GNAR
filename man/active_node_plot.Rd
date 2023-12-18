\name{active_node_plot}
\alias{active_node_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces an active node matrix heat-map.
%%  ~~function to do ... ~~
}
\description{Produces an active node matrix heat-map, which compares the local impact each node has on all the other ones (i.e., regressing \eqn{j} on \eqn{i}) once a model order has been chosen. The local relevance indes is 
\eqn{\mathrm{local} (i, j) := \bigg ( w_{ij} \sum_{k = 1}^{p} |\hat{\beta}_{kr}| \bigg ) \bigg \{ \sum_{l \in \mathcal{N} (i)} \sum_{r = 1}^{r^*} \sum_{k = 1}^{p}  w_{il} |\hat{\beta}_{kr}| \bigg) \bigg \}^{-1},}
which is closer to one the more relevant \eqn{j} is when forecasting \eqn{i}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  active_node_plot(vts, network, max_lag, r_stages)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vts}{Vector time series under study.}
  \item{network}{GNAR network object, which is the underlying network for the time series under study.}
  \item{max_lag}{Maximum lag of the fitted GNAR model - i.e., \eqn{\mathrm{GNAR}(p, [s_1, \dots, s_p]).}}
  \item{r_stages}{Neighbourhood regression oreder of the fitted GNAR model - i.e., \eqn{(s_1, \dots, s_p)}.}
%%     ~~Describe \code{x} here~~
}
\value{Produces the local influence matrix heat-map for a specific model order. Does not return any values.
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
# Produces an active node heat-map matrix from a stationary GNAR(2, [2, 1]) simulation.
#
gnar_simulation <- GNARsim(n = 100, net=fiveNet,
	alphaParams = list(rep(0.25, 5), rep(0.12, 5)), 
        betaParams = list(c(0.25, 0.13), c(0.20)), sigma=1)
#
# Active node plot
#
active_node_plot(gnar_simulation, fiveNet, 2, c(2, 1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
