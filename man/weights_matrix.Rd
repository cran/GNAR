\name{weights_matrix}
\alias{weights_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes the weights matrix corresponding to the GNAR network object linked to the vector time series.
%%  ~~function to do ... ~~
}
\description{Computes the weights matrix with normalised weights (i.e., add up to one) for the network time series with underlying network provided by the user. If the network is unweighted, then each r-stage neighbour is considered equally relevant, i.e., \eqn{w_{ij}} = \eqn{\{ \mathcal{N}_r (i)\}^{-1} \mathbb{I} (d(i, j) = r) }, where \eqn{\mathbb{I}} is the indicator function and the distance is the shortest path in the underlying network.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
weights_matrix(network, max_r_stage)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{Network linked to the vector time series under study, must be a GNARnet object.}
  \item{max_r_stage}{Longest shortest path for which weights are non-zero. If not specified, then its set equal to the upper bound, which is the longest shortest path in the underlying network.}
  
%%     ~~Describe \code{x} here~~
}
\value{Weight matrix \eqn{\mathbf{W}}, each entry is the weight \eqn{w_{ij}} between a pair of nodes. The matrix is not symmetric, and each row adds up to one when considering r-stage neighbours for a particular r.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
# Weights matrix linked to the mechanical ventilation beds time series. 
# This network has a longest shortest path equal to six.
# 
#data(fiveNet)
W_norm = weights_matrix(fiveNet, 6)
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
