\name{summary.GNARfit}
\alias{summary.GNARfit}
\title{Returns model summary for a GNAR model fit}
\usage{
\method{summary}{GNARfit}(object, ...)
}
\description{
Returns the summary of a GNARfit object, including BIC.
}
\arguments{
\item{object}{output of a \link[GNAR]{GNARfit} call.}
\item{...}{additional arguments, unused here.}
}

\details{
The output is the summary of the fit using \link[stats]{summary.lm}, and BIC calculated using \link[GNAR]{BIC.GNARfit}.
}
\value{
\code{summary.GNARfit} prints the model summary and the value of the BIC.
}
\examples{
#summary for the GNAR(2,[1,1]) model using GNARfit on fiveVTS
summary(GNARfit())
}
