\name{seed.nos}
\alias{seed.nos}
\title{Vector of seed numbers}
\usage{
seed.nos
}
\description{
Seed numbers for reproducible random graphs.
}
\format{
\code{seed.nos} is a vector of length 10,000 containing integers.
}
\examples{
data(gdp)
g <- seedToNet(seed.nos[1], nnodes=35, graph.prob=0.15)
plot(g, vertex.label=colnames(gdpVTS), vertex.size=0)
}
