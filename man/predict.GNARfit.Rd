\name{predict.GNARfit}
\alias{predict.GNARfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction of a GNARfit object
}
\description{
Predicts future observations from a GNARfit object, based on the fitted GNAR model.
}
\usage{
\method{predict}{GNARfit}(object, n.ahead=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the output of a \link[GNAR]{GNARfit} call
}
  \item{n.ahead}{the time length of the predictions}
  \item{...}{further arguments passed to the \link[GNAR]{simulate.GNARfit} function, such as \code{seed} }
}
\details{
S3 method for class "GNARfit". This function calls \link[GNAR]{simulate.GNARfit}.
}
\value{
A multivariate time series of dimension \code{n.ahead} x nnodes.
}

\examples{
#simulate 5 future observations from fiveVTS
predict(GNARfit(), n.ahead=5)
}
