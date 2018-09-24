\name{NofNeighbours}
\alias{NofNeighbours}
\title{Calculates stage-neighbours of a network}
\usage{
NofNeighbours(node=1, stage=2, net=GNAR::fiveNet)
}
\description{
Calculates neighbour sets of a particular node in the network and their distances.
}
\arguments{
\item{node}{is an integer specifying which node to calculate the neighbours of.}
\item{stage}{is an integer specifying the maximum neighbour-stage to calculate to.}
\item{net}{a \code{GNARnet} object with edge list and distance list.}
}
\details{
Note that the distances are calculated as the sum along the shortest path; do not use this with a weights (rather than distance) list. Stage-\code{r} neighbours of node \code{i} are denoted \eqn{N^(r)(i)}{N^(r) (i)}, and are nodes that are \code{r} edges (but no fewer) away from \code{i}. Hence stage-1 neighbours are the immediate neighbours, stage-2 neighbours are the neighbours of neighbours and so on.
}
\value{
\item{edges}{is a list of length \code{stage}, where \code{edges[[i]]} is a vector containing the nodes that are stage-\code{i} neighbours of \code{node}.}
\item{dist}{is a list of length \code{stage}, where \code{dist[[i]]} is a vector containing the distances from \code{node} to its stage-\code{i} neighbours, with ordering as in \code{edges[[i]]}.}
}
\examples{
#First and second stage neighbours of node 1 in fiveNet
data(fiveNode)
NofNeighbours()
}
