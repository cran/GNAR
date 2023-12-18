\name{BIC.GNARfit}
\alias{BIC.GNARfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Information Criterion for GNAR models
}
\description{
Function calculating BIC for \code{GNARfit} models.
}
\usage{
\method{BIC}{GNARfit}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{GNARfit} object, output from a \link[GNAR]{GNARfit} call.}
  \item{...}{additional arguments, not used here.}
}
\details{
Smaller BIC values correspond to better fit.
}
\value{
A numeric value corresponding to the \code{BIC}.  Note that this is the ``time-normalised'' value of the AIC for the GNAR model, and also removes any proportionality constants in the calculation.
}
\examples{
#BIC for two different GNAR fits for fiveNet data
#GNAR(2,[1,1])
BIC(GNARfit())
#GNAR(2,[1,0])
BIC(GNARfit(betaOrder=c(1,0)))
}

