\name{corbit_plot}
\alias{corbit_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Corbit (correlation-orbit) plot, which aids model selection by visualising network autocorrelation and partial network autocorrelation.
%%  ~~function to do ... ~~
}
\description{Plots the GNAR network autcorrelation funciton for a choice of maximum lag and maximum r-stage depth in the network. Using the \code{\link[GNAR]{nacf}} function for network autocorrelation and \code{\link[GNAR]{pnacf}} for partial network autocorrelation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
corbit_plot(vts, net, max_lag, max_stage, weight_matrix, 
viridis_color_option="viridis", size_option="absolute_val", 
partial="no", r_corbit="no", line_segment="no", rectangular_plot="no")
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
  \item{size_option}{Point size scale for the Corbit plot. The default is the absolute value of the network autocorrelation function (i.e., \code{|nacf(h, r)|} or  \code{|pnacf(h, r)|)}. Alternate option is the coefficient of determination coming from a global-\eqn{\alpha} model with fixed lag and stage.}
  \item{partial}{Option for selecting between computing the network autocorrelation function or the partial network autocorrelation funciton. Default choice is network autocorrelation (i.e., partial="no"), change argument to "yes" for computing the partial network autocorrelation funciton.}
  \item{r_corbit}{Choice for distinguishing between Corbit and R-Corbit plots, default is set to Corbit (inner function call). For producing R-Corbit plots one should use \code{\link[GNAR]{r_corbit_plot}}. }
   \item{line_segment}{Default is set to no (i.e., \code{line_segment = "no"}). Changing the option to yes produces an alternative Corbit plot with a line segment between the first and last lags. }
  \item{rectangular_plot}{Option for producing alternative rectangular plots. Default is set to no (i.e., \code{rectangular_plot = "no"}). Changing it to yes (i.e., \code{rectangular_plot = "yes"}) produces a rectangular plot with (P)NACF values on the y-axis and lags on the x-axis, r-stages are distinguished by colour. Changing it to square (i.e., \code{rectangular_plot = "square"}) produces a square plot with r-stages on the y-axis and lags on the x-axis.}

}
\details{This function computes network autocorrelation or partial network autocorrelation function values for a specific choice of maximum lag and r-stage depth, and produces the corresponing Corbit plot. Each point in a Corbit plot corresponds to the (P)NACF value, i.e., \code{(p)nacf(h, r)}, at a h-lag and r-stage pair. The ring number starting from the inside corresponds to r-stage depth (path length), and the numbers on the outside ring indicate time lag. The colour scale is based on the overall network autocorrelation values (i.e., the colour is set to highlight strong correlations). The dot at the centre has \code{(p)nacf(h, r)=0} and the smallest size for aiding comparison. Please see the RMarkdown document mentioned in \code{GNAR} references for plot examples.
%%  ~~ If necessary, more details than the description above ~~
}

\value{Produces the specified, i.e, NACF or PNACF for a choice of lag and \eqn{r}-stage depth, \eqn{(h, r)}, Corbit plot. Does not print (P)NACF values, these are stored as an invisble data frame (matrix), and can be accessed by printing or calling the object produced by the \code{corbit_plot} call.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Nason, G.P., Salnikov, D. and Cortina-Borja, M. (2023)
New tools for network time series with an application to COVID-19
hospitalisations. \url{https://arxiv.org/abs/2312.00530}
}
\author{Daniel Salnikov and Guy Nason.
%%  ~~who you are~~
}

\examples{
#
# Simulate 100 observations from a stationary GNAR(2, [2, 1]), where 
# fiveNet is the underlying network.
# 
gnar_simulation <- GNARsim(n = 100, net=fiveNet, alphaParams = list(rep(0.25, 5), rep(0.12, 5)), 
        betaParams = list(c(0.25, 0.13), c(0.20)), sigma=1)
# We produce the corresponding Corbit plots.
corbit_plot(gnar_simulation, fiveNet, 20, 3)
corbit_plot(gnar_simulation, fiveNet, 20, 3, partial = "yes")

# If the network object comes with its own weights, then these can be added by including 
# the option weight_matrix in the corbit call.
# corbit_plot(vts, net, max_lag, max_stage, weight_matrix = object_weights_matrix)
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
