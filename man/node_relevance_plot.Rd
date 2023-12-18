\name{node_relevance_plot}
\alias{node_relevance_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces a node relevance plot, which compares the impact each node has on the network autocorrelation once a model order has been chosen.
%%  ~~function to do ... ~~
}
\description{Produces a node relevance plot based on the node relevance index \eqn{\mathrm{globindex}(X_{i, t}) := \bigg (\sum_{j = 1}^{d} [\mathbf{W} \odot \mathbf{S}]_{j i} \bigg )
   \bigg \{ \underset{l \in \mathcal{K}}{\max} \bigg ( \sum_{j = 1}^{d} [\mathbf{W} \odot \mathbf{S}]_{j l} \bigg ) \bigg \}^{-1}, } which computes the ratio between nodes \eqn{i} column sum for nodes in neighbourhood regressions. Nodes are ordered according to the relative contribution eahc has to the autocovariance. The nodes are ordered in ascending order.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  node_relevance_plot(network, r_star, node_names, node_label_size = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{GNAR network object, which is the underlying network for the time series under study.}
  \item{r_star}{Maximum active r-stage depth for neighbourhood regression.}
  \item{node_names}{Names corresponding to each, this makes identifying nodes in the plot easier. If this argument is NULL, then the plot links to each node a number.}
  \item{node_label_size}{Text size when producing the plot. Default is 2, however, depending on the number of nodes it might be necessary to adjust the size.}
}
\value{
Data Frame consisting of two variable, the node name and the node relevance value.
}
\author{Daniel Salnikov and Guy Nason.
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
# Produces a node relevance plot with respect to a stationary GNAR process 
# with underlying network fiveNet
#
# GNAR simulation
gnar_simulation <- GNARsim(n = 100, net=fiveNet, alphaParams = list(rep(0.25, 5), rep(0.12, 5)), 
        betaParams = list(c(0.25, 0.13), c(0.20)), sigma=1)
# Node relevance plot without names
node_relevance_plot(network = fiveNet, r_star = 2, node_label_size = 10)
#
# Node relevance plot with names
#
node_relevance_plot(network = fiveNet, r_star = 2, node_names = c("A", "B", "C", "D", "E"), 
	node_label_size = 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ nacf }
\keyword{ relevance }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ GNAR }
% \concept{ ~cpt2 }
% Use only one concept per line.
