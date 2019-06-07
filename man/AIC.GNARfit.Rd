\name{AIC.GNARfit}
\alias{AIC.GNARfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Akaike's Information Criterion for GNAR models
}
\description{
Function calculating AIC for \code{GNARfit} models.
}
\usage{
\method{AIC}{GNARfit}(object, ..., k=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{GNARfit} object, output from a \link[GNAR]{GNARfit} call.}
  \item{...}{additional arguments, not used here.}
  \item{k}{the penalty for the criterion, the default \code{k = 2} is the standard AIC.}
}
\details{
Smaller AIC values correspond to better fit.
}
\value{
A numeric value corresponding to the \code{AIC} (or other criterion if \code{k} is set to something other than \code{2}).
}
\examples{
#AIC for two different GNAR fits for fiveNet data
#GNAR(2,[1,1])
AIC(GNARfit())
#GNAR(2,[1,0])
AIC(GNARfit(betaOrder=c(1,0)))
}

