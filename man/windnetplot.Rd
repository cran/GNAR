\name{windnetplot}
\alias{windnetplot}
\title{Produce bespoke plot of the wind data network}
\usage{
windnetplot()
}
\description{
Plots the wind speed data network with distance information.
}
\arguments{
	None.
}
\details{The wind speed data is to be found in the \code{\link{vswind}}
	data set. This function contains commands, using 
	functionality from the \code{wordcloud} package, to plot the
	network, with node names and edges. Distances between nodes are
	plotted next to the edges.
}
\seealso{\code{\link{vswind}}}
\examples{
data(vswind)
\dontrun{windplotnet()}
}
