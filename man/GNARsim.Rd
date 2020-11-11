\name{GNARsim}
\alias{GNARsim}
\title{Simulates a GNAR process}
\usage{
GNARsim(n=200, net=GNAR::fiveNet, alphaParams=list(c(rep(0.2,5))),
 betaParams=list(c(0.5)), sigma=1, tvnets=NULL, netsstart=NULL)
}
\description{
Simulates a GNAR process with Normally distributed innovations.
}
\arguments{
\item{n}{time length of simulation.}
\item{net}{network used for the GNAR simulation.}
\item{alphaParams}{a list containing vectors of auto-regression parameters for each time-lag.}
\item{betaParams}{a list of equal length as \code{alphaParams} containing the network-regression parameters for each time-lag.}
\item{sigma}{the standard deviation for the innovations.}
\item{tvnets}{Only NULL is currently supported.}
\item{netsstart}{Only NULL is currently supported.}
}
\details{
Parameter lists should not be NULL, set unused parameters to be zero. See \link[GNAR]{GNARfit} for model description.
}
\value{
\code{GNARsim} returns the multivariate time series as a \link[stats]{ts} object, with \code{n} rows and a column for each of the nodes in the network.
}
\references{Knight, M.I., Nunes, M.A. and Nason, G.P. Modelling, detrending and decorrelation of network time series.
	\href{https://arxiv.org/abs/1603.03221}{arXiv preprint}.\cr\cr
Knight, M.I., Leeming, K., Nason, G.P. and Nunes, M. A. (2020) Generalised Network Autoregressive Processes and the GNAR package. \emph{Journal of Statistical Software}, \bold{96} (5), 1--36.
}
\examples{
#Simulate a GNAR(1,[1]) process with the fiveNet network
GNARsim()
}
