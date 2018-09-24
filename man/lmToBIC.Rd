\name{lmToBIC}
\alias{lmToBIC}
\title{Calculates the BIC for a GNAR model fit}
\usage{
lmToBIC(GNARobj.in=GNARpredict())
}
\description{
Returns the value of the BIC for a GNARfit or GNARpredict object.
}
\arguments{
\item{GNARobj.in}{output of \code{GNARfit} or \code{GNARpredict}.}
}
\details{
For a GNAR(p,[S]) model, the BIC is calculated as
\deqn{BIC(p,S) = \ln | \Sigma_{p,S} | + T^{-1} M\ln(T)}{BIC(p,S) = log |Sigma_{p,S}| + 1/T M log(T)}
where \eqn{\Sigma_{p,S}=T^{-1}U'U}{Sigma_{p,S}=1/T U'U}, \eqn{U} is the matrix of residuals, and \eqn{M=np + C\sum_{j=1}^p S_j}{M=np+C(S_1 + ... + S_p)} when \code{globalalpha=FALSE},  \eqn{M=p + C\sum_{j=1}^p S_j}{M=p+C(S_1 + ... + S_p)} when \code{globalalpha=TRUE}, and \code{C} is the number of factors.
}
\value{
\code{lmToBIC} returns the value of the BIC for the input model fit.
}
\examples{
#BIC for the GNAR(2,[1,1]) model using GNARpredict on fiveVTS
data(fiveNode)
lmToBIC()
}
