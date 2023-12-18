\name{wagner_plot}
\alias{wagner_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces a Wagner plot for the specified choice of covariates and/or time slices.
}
\description{Produces a Wagner plot for comparing the netowrk autocorrelation and partial network autocorrerlation function values for a choice of maximum lag and maximum r-stage depth. The plot consists of rings linked to a time-slice or covariate, and indicates if the network autocorrelation differs for different time-slices and/or covariates. It is a visual for studying whether or not the network time series acts differently depending on lag or covariate. Starting from one and continuing to the outside, each ring corresponds to said choice of r-stage depth, the numbers on the outside ring are time-lags, and each dot corresponds to a specific time-slice or covariate.
}
\usage{
wagner_plot(vts_frames, network_list, max_lag, max_stage, weight_matrices, 
frame_names, same_net="no", viridis_color_option="viridis", size_option="absolute_val", 
partial="no", wagner="yes")
}
\arguments{
  \item{vts_frames}{List containing the vector time series linked to each of the covariate and/or time-slice, which the Wagner plot compares.}
  \item{network_list}{List of network objects for which the Wagner plot compares network autocorrelation or partial network autocorrelation.}
  \item{max_lag}{Maximum lag for the Wagner plot.}
  \item{max_stage}{Maximum r-stage depth for the Wagner plot (i.e., the number of rings in the Wagner plot).}
  \item{weight_matrices}{List of weight matrices, each weight matrix corresponds to a particular choice of time-slice or covariate. If all the time-slices have the same weight matrix, then the argument is a list with said list.}
  \item{frame_names}{Indicates the name of each time-slice or covariate time series, order should be the same as in the weight matrices and vector time series lists.}
  \item{same_net}{Indicates whether or not all time-slices or covariates share the same weight matrix. Default choice is no, if the time-slices or covariates share the same weight matrix, then this argument should be set to "yes" (i.e., same_net = "yes").}
  \item{viridis_color_option}{Colour scale for the Corbit plot. The default option is \code{viridis}, each option is colout blind friendly; see \code{viridis} package.}
  \item{size_option}{Point size scale for the Corbit plot. The default is the absolute value of the network autocorrelation function (i.e., \code{|nacf(h, r)|} or  \code{|pnacf(h, r)|}). Alternate option is the coefficient of determination coming from a global-\eqn{\alpha} model with fixed lag and stage.}
   \item{partial}{Option for selecting between computing the network autocorrelation function or the partial network autocorrelation funciton. Default choice is network autocorrelation (i.e., partial="no"), change argument to "yes" for computing the partial network autocorrelation funciton.}
    \item{wagner}{Choice for distinguishing between Corbit and Wagner plots, default is set to Corbit (inner function call). For producing Wagner plots one should use \code{\link{wagner_plot}}. }
}
\details{Wagner plots compare the network autocorrelation function (NACF) and partial network autocorrelation function (PNACF) values of different time-slices and/or covariate weights. Wagner plots are read in the same manner as Corbit plots \code{\link{corbit_plot}}, and include a legend on the right-hand side for distinguishing between covariates and/or time-slices. The point at the centre is the mean value of the NACF or PNACF values arising from the time-slices and/or covariate data splits. Essentially, if \eqn{c \in \{1, \dots, C\}}, where \eqn{C \in \mathbb{N}} is the number of covariates or time-slices, then the value at the centre is
\eqn{\mathrm{nacf}(h, r) = C^{-1} \sum_{c = 1}^{C} \mathrm{nacf}_c (h, r),} where \eqn{\mathrm{nacf}_c(h, r)} is the NACF value corresponding to the covariate/time-slice \eqn{c}. The number of covariates and time-slices \eqn{C} must be equal to the length of the lists used for producing the Wagner plot.
}
\value{Produces the specified Wagner plot, does not return the network autocorrelaiton values.
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
# Produces a Wagner plot, which compares three stationary GNAR simulations, 
# where the underlying network is fiveNet.
#
# Compute the weight matrix
W = weights_matrix(fiveNet)
#
# Simulate three stationary GNAR processe
sim1 <-  GNARsim(n = 100, net=fiveNet, alphaParams = list(c(0.1, 0.12, 0.16, 0.075, 0.21),
                                                          c(0.12, 0.14, 0.15, 0.6, 0.22)), 
                 betaParams = list(c(0.1, 0.16), c(0.11, 0.14)))

sim2 <-  GNARsim(n = 100, net=fiveNet, alphaParams = list(rep(.25, 5)), 
                 betaParams = list(c(0.1, 0.16)))

sim3 <- GNARsim(n = 100, net=fiveNet, alphaParams = list(rep(.25, 5), rep(0.13, 5)), 
                betaParams = list(c(0.1, 0.16), c(0.11)))
                
# Produce NACF Wagner plot with the same network and weights matrix
wagner_plot(list(sim1, sim2, sim3), list(fiveNet), 10, 3, list(W), 
		c("sim1", "sim2", "sim3"), same_net = "yes")
# 
# Produce PNACF Wagner with different networks and weight matrices
wagner_plot(list(sim1, sim2, sim3), list(fiveNet, fiveNet, fiveNet), 10, 3, list(W, W, W), 
		c("sim1", "sim2", "sim3"), same_net = "no", partial = "yes")
}

}


\keyword{ nacf }
\keyword{ pnacf }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ GNAR }
\concept{ Wagner }
\concept{ Corbit }
% Use only one concept per line.
