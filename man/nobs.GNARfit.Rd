\name{nobs.GNARfit}
\alias{nobs.GNARfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to return the number of observations input to GNARfit objects
}
\description{
\code{nobs} returns the number of obervations (T) of the input multivariate time series in the \code{GNARfit} function.
}
\usage{
\method{nobs}{GNARfit}(object,...)
}
\arguments{
  \item{object}{
the output of a GNARfit or GNARpredict call
}
\item{...}{additional arguments, unused here.}
}
\details{
S3 method for class "GNARfit".
}
\value{
An integer specifying the number of rows in the input \code{vts} to the \link[GNAR]{GNARfit} function.
}

\examples{
#observations of example fiveVTS
nobs(GNARfit())
#check this is the same as number of rows in fiveVTS
all.equal(nobs(GNARfit()), nrow(fiveVTS))
}
