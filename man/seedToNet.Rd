\name{seedToNet}
\alias{seedToNet}
\title{Produces a random network from a seed value}
\usage{
seedToNet(seed.no, nnodes=34, graph.prob=0.5)
}
\description{
Produces a reproducible undirected Erdos-Reyni random network using a particular seed value.
}
\arguments{
\item{seed.no}{a valid number to set the seed to.}
\item{nnodes}{the number of nodes in the produced network.}
\item{graph.prob}{the probability that each pair of nodes is connected.}
}
\details{
\code{graph.prob} effectively controls the sparsity of the network. All distances are set to 1.
}
\value{
A GNARnet object.
}
\examples{
#Generate the random graph from seed 10, with 5 nodes and connection prob 0.5
seedToNet(10,nnodes=5,graph.prob=0.5)
}
