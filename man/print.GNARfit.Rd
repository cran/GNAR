\name{print.GNARfit}
\alias{print.GNARfit}
\title{Function to print the model and coefficients of GNARfit objects}
\usage{
\method{print}{GNARfit}(x,...)
}
\description{
\code{print.GNARfit} prints model, call, and coefficients of a GNARfit object.
}
\arguments{
\item{x}{the output of a \link[GNAR]{GNARfit} call}
\item{...}{additional arguments, unused here.}
}
\details{
S3 method for class "GNARfit".
}
\examples{
#print the information of the fiveNode GNAR fit
print(GNARfit())
}
