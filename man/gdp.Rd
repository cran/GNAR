\name{gdp}
\alias{gdp}
\docType{data}
\title{
GDP example data.
}
\description{
Data and seed values for GDP example.
}
\usage{data("gdp")}
\format{ This dataset contains two R objects:\cr
\code{gdpVTS} is a \link[stats]{ts} object with a matrix of 52 rows (t=52) and 35 columns (n=35)
\code{seed.nos} is a vector containing integer values.
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
data(gdp)
image(gdpVTS)
plot(seedToNet(seed.nos[1], nnodes=35, graph.prob=0.15))
}
\keyword{datasets}
