\name{simulate.GNARfit}
\alias{simulate.GNARfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate from a GNARfit object
}
\description{
Simulates from a GNARfit object, either creating a new series or future observations of the original series based upon the fitted GNAR model.
}
\usage{
\method{simulate}{GNARfit}(object, nsim=object$frbic$time.in, seed=NULL,
                             future=TRUE, set.noise=NULL, allcoefs=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the output of a \link[GNAR]{GNARfit} call
}
  \item{nsim}{the time length of the simulations}
  \item{seed}{either \code{NULL}, or a value to set the seed to}
  \item{future}{whether the simulations follow on from the original time series (\code{TRUE}), or if \code{FALSE} the simulations are a new series.}
  \item{set.noise}{a value to set the standard deviation of the noise to, or if \code{NULL}, the estimated standard deviation from the input series will be used.}
  \item{allcoefs}{if \code{TRUE}, all fitted coefficients will be used, if \code{FALSE} only the significant (p-val < 0.05) coefficients will be used.}
  \item{...}{additional arguments, unused here.}
}
\details{
S3 method for class "GNARfit".
}
\value{
A multivariate time series of dimension \code{nsim} x nnodes.
}

\examples{
#simulate 5 future observations from fiveVTS
simulate(GNARfit(), nsim=5)
}
