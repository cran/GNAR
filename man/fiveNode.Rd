\name{fiveNode}
\alias{fiveNode}
\alias{fiveVTS}
\alias{fiveNet}
\docType{data}
\title{
Example Network Time Series
}
\description{
A multivariate time series \code{fiveVTS} and corresponding network \code{fiveNet}.
}
\usage{data("fiveNode")}
\format{ This dataset contains two R objects:\cr
\code{fiveVTS} is a \link[stats]{ts} object with a matrix of 200 rows (t=200) and 5 columns (n=5)
\code{fiveNet} is a GNARnet object containing \code{$edges} and \code{$dist}. \cr
\code{edges} is a list of length five, with \code{edges[[i]]} containing the vertices that node i is connected to. \cr
\code{dist} is a list of length five, with \code{dist[[i]]} containing the length of the vertices that node i is connected to.

}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
plot(fiveNet)
image(fiveVTS)
}
\keyword{datasets}
