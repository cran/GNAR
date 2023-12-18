\name{get_k_stages_adjacency_tensor}
\alias{get_k_stages_adjacency_tensor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes a list of r-stage adjacency matrices.
%%  ~~function to do ... ~~
}
\description{Computes a list of r-stage adjacency matrices, each matrix in the list inidicates whether or not nodes \eqn{i} and \eqn{j} are r-stage neighbours in the underlying network. Essentially \eqn{[\mathbf{S}_r]_{ij} = 1 } if and only if \eqn{d(i, j) = r}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  get_k_stages_adjacency_tensor(St_1, r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{St_1}{One-stage adjacency matrix (i.e., the adjacency matrix of the underlying network).}
  \item{r}{Maximum r-stage for which one wishes to compute the r-stage adjacency matrix.}
%%     ~~Describe \code{x} here~~
}
\value{List containing the adjacency matrices in ascending order.
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{\eqn{ \{ \mathbf{S}_q \}_{q=1}^{q = r} }}{Each entry is the r-stage adjacency matrix at depth \eqn{r}.}
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
# Produce the r-stage adjacency tensors for the fiveNet network.
#
get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(fiveNet)), 3)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ rstage }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
