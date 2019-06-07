\name{vcov.GNARfit}
\alias{vcov.GNARfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate variance-covariance matrix for a itted GNARfit object
}
\description{
Returns the varaince-covariance matrix of the parameters of a GNARfit object.
}
\usage{
\method{vcov}{GNARfit}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a GNARfit object, the output from a \link[GNAR]{GNARfit} call.
}
  \item{...}{further arguments passed to \link[stats]{vcov}.}
}
\details{
S3 method for class "GNARfit".
}
\value{
A matrix of estimated covariances between the parameter estimates, this is calculated using \link[stats]{vcov} for \link[stats]{lm} objects.
}

\examples{
#covariance matrix of fiveNode fit
vcov(GNARfit())
}
