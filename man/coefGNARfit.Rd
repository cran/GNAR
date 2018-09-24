\name{coef.GNARfit}
\alias{coef.GNARfit}
\title{Function to return coefficients of GNARfit objects}
\usage{
\method{coef}{GNARfit}(object,...)
}
\description{
\code{coef.GNARfit} returns the vector of coefficients from a GNARfit object.
}
\arguments{
\item{object}{the output of a GNARfit or GNARpredict call}
\item{...}{additional arguments, unused here.}
}
\details{
S3 method for class "GNARfit".
}
\value{
\code{coef.GNARfit} returns a vector of coefficient values.
}
\examples{
#get the coefficients of the fiveNode data GNAR fit
data(fiveNode)
coef(GNARfit())
}
