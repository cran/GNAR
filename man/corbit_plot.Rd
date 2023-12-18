\name{corbit_plot}
\alias{corbit_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Corbit (correlation-orbit) plot, which aids model selection by visualising network autocorrelation and partial network autocorrelation.
%%  ~~function to do ... ~~
}
\description{Plots the GNAR network autcorrelation funciton for a choice of maximum lag and maximum r-stage depth in the network. Using the \code{\link{nacf}} function for network autocorrelation and \code{\link{pnacf}} for partial network autocorrelation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
corbit_plot(vts, net, max_lag, max_stage, weight_matrix, 
    viridis_color_option="viridis", size_option="absolute_val", 
    partial="no", wagner="no")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vts}{Vector time series observations for which one wishes to plot the network autocorrelation or partial network autocorrelation.
%%     ~~Describe \code{x} here~~
  }
  \item{net}{GNAR network object linked to the time series under study.}
  \item{max_lag}{Maximum lag the Corbit plot produces (i.e., number of time-steps considered for the network autocorrelaiton.)}
  \item{max_stage}{Maximum r-stage depth considered for the Corbit plot (i.e., the number of rings in the plot). Corresponds to the length of paths in the underlying network.}
  \item{weight_matrix}{A matrix which entries correspond to the weights between nodes. If this term is NULL, then this argument is equal weights between r-stage neighbours.}
  \item{viridis_color_option}{Colour scale for the Corbit plot. The default option is \code{viridis}, each option is colout blind friendly; see \code{viridis} package.}
  \item{size_option}{Point size scale for the Corbit plot. The default is the absolute value of the network autocorrelation function (i.e., \code{|nacf(h, r)|} or  \code{|pnacf(h, r)|}). Alternate option is the coefficient of determination coming from a global-\eqn{\alpha} model with fixed lag and stage.}
  \item{partial}{Option for selecting between computing the network autocorrelation function or the partial network autocorrelation funciton. Default choice is network autocorrelation (i.e., partial="no"), change argument to "yes" for computing the partial network autocorrelation funciton.}
  \item{wagner}{Choice for distinguishing between Corbit and Wagner plots, default is set to Corbit (inner function call). For producing Wagner plots one should use \code{\link{wagner_plot}}. }

}
\details{Function calculates the network autocorrelation or the partial network autocorrelation values for a specific choice of maximum lag and r-stage depth, and produces the corresponing Corbit plot. Each point in the Corbit plot corresponds to the network autocorrelation \code{nacf(h, r)} at a h-lag and r-stage pair. The ring number starting from the inside corresponds to r-stage depth (path length), and the numbers on the outside ring indicate the time lag. The colour scale is based on the overall network autocorrelation values (i.e., the colour is set to highlight strong correlations).
%%  ~~ If necessary, more details than the description above ~~
}
\value{Produces the specified Corbit plot, does not return the network autocorrelaiton values.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Daniel Salnikov and Guy Nason.
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\references{
Nason, G.P., Salnikov, D. and Cortina-Borja, M. (2023)
New tools for network time series with an application to COVID-19
hospitalisations. \url{https://arxiv.org/abs/2312.00530}
}

\examples{

\dontrun{
#
# Simulate 100 observations from a stationary GNAR(2, [2, 1]), where fiveNet is 
# the underlying network.
# 
gnar_simulation <- GNARsim(n = 100, net=fiveNet, alphaParams = list(rep(0.25, 5), rep(0.12, 5)), 
        betaParams = list(c(0.25, 0.13), c(0.20)), sigma=1)
# We produce the corresponding Corbit plots.
corbit_plot(gnar_simulation, fiveNet, 20, 3)
corbit_plot(gnar_simulation, fiveNet, 20, 3, partial = "yes")

# If the network object comes with its own weights, then these can be added by including the 
# option weigh_matrix in the corbit call.
# corbit_plot(vts, net, max_lag, max_stage, weight_matrix = object_weights_matrix)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ corbit }
\keyword{ nacf }
\keyword{ pnacf }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
