\name{is.GNARfit}
\alias{is.GNARfit}
\title{Function to check GNARfit objects}
\usage{
is.GNARfit(x)
}
\description{
\code{is.GNARfit} returns either TRUE or FALSE according to a series of GNARfit checks.
}
\arguments{
\item{x}{the object to be tested}
}
\details{
The \code{is.GNARfit} function checks whether the object passes a series of tests that correspond to it being the output of \link[GNAR]{GNARfit} or \link[GNAR]{GNARpredict}:
\itemize{
  \item Is it a list containing \code{$mod} and \code{$frbic}
  \item Does it contain either \code{$y} and \code{$dd} or $\code{ys} and \code{$ds}
  \item Is \code{$mod} a \link[stats]{lm} object
  \item Does \code{$frbic} have the components to calculate the BIC with \link[GNAR]{lmToBIC}
}
}
\value{
\code{is.GNARfit} returns \code{TRUE} or \code{FALSE} corresponding to passing the above tests.
}
\examples{
#check that the example fit meets the criteria above
data(fiveNode)
is.GNARfit(GNARfit())
#also check the predict function
is.GNARfit(GNARpredict())
}
