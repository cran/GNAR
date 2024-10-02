\name{r_corbit_plot}
\alias{r_corbit_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces a R-Corbit plot for the specified choice of covariate-levels and/or time-slices.
}
\description{Produces a R-Corbit plot for comparing the network autocorrelation (NACF) and partial network autocorrelation function (PNACF) values for a choice of maximum lag and maximum \eqn{r}-stage depth. Starting from the first and continuing to the outermost ring, each ring corresponds to said choice of \eqn{r}-stage depth. The numbers on the outermost ring are time-lags, and each dot corresponds to a specific time-slice or covariate-level.
}
\usage{
r_corbit_plot(vts_frames, network_list, max_lag, max_stage, weight_matrices, 
frame_names, same_net="no", viridis_color_option="viridis", size_option="absolute_val", 
partial="no", r_corbit="yes")
}
\arguments{
  \item{vts_frames}{List containing the vector time series linked to each of the covariate-levels and/or time-slices, which the R-Corbit plot compares.}
  \item{network_list}{List of network objects for which the R-Corbit plot compares network autocorrelation or partial network autocorrelation.}
  \item{max_lag}{Maximum lag for the R-Corbit plot.}
  \item{max_stage}{Maximum \eqn{r}-stage depth for the R-Corbit plot (i.e., the number of rings in the R-Corbit plot).}
  \item{weight_matrices}{List of weigth matrices, each weight matrix corresponds to a particular choice of time-slice or covariate-level. If all the time-slices have the same weight matrix, then the argument is a list, where all the entries are equal to the unique weight matrix.}
  \item{frame_names}{Indicates the name of each time-slice or covariate-level time series. Order should be the same as in the weight matrices and vector time series lists.}
  \item{same_net}{Indicates whether or not all time-slices or covariate-levels share the same weight matrix. Default choice is no, if the time-slices or covariate-levels share the same weight matrix, then this argument should be set to "yes" (i.e., same_net = "yes").}
  \item{viridis_color_option}{Colour scale for the R-Corbit plot. The default option is \code{viridis}, each option is colout blind friendly.}
  \item{size_option}{Point size scale for the R-Corbit plot. Default is the absolute value of the network autocorrelation function (i.e., \code{|nacf(h, r)|} or  \code{|pnacf(h, r)|)}. Alternate option is the coefficient of determination coming from a global-\eqn{\alpha} model constrained to a specific lag and stage pair.}
   \item{partial}{Option for selecting between computing the network autocorrelation function or the partial network autocorrelation function. Default choice is network autocorrelation (i.e., partial="no"). Change argument to "yes" for computing the partial network autocorrelation function (PNACF).}
    \item{r_corbit}{Choice for distinguishing between Corbit and R-Corbit plots, default is set to Corbit (inner function call). For producing R-Corbit plots one should use \code{\link[GNAR]{corbit_plot}}. }
}
\details{R-Corbit plots compare the network autocorrelation function (NACF) and partial network autocorrelation function (PNACF) values for a choice of different time-slices and/or covariate-levels. R-Corbit plots are read in the same manner as Corbit plots \code{\link{corbit_plot}}, and include a legend on the right-hand side for distinguishing between covariate-levels and/or time-slices. The point at the centre is the mean value of the NACF or PNACF values arising from the time-slices and/or covariate-levels data splits. Essentially, if \eqn{c \in \{1, \dots, C\}}, where \eqn{C \in \mathbb{N}} is the number of covariate-levels or time-slices, then the value at the centre is \eqn{\mathrm{(p)nacf}(h, r) = C^{-1} \sum_{c = 1}^{C} \mathrm{(p)nacf}_c (h, r),} where \eqn{\mathrm{(p)nacf}_c(h, r)} is the (P)NACF value corresponding to the covariate-level/time-slice \eqn{c}. The number of covariate-levels and time-slices \eqn{C} must be equal to the length of the lists used for producing the R-Corbit plot.
}
\value{Produces the specified, i.e., NACF or PNACF, values for a choice of lag and \eqn{r}-stage depth, \eqn{(h, r)}, R-Corbit plot. Does not print (P)NACF values, these are stored as invisble data frames (matrices), and can be accessed by printing or calling the object produced by the \code{r_corbit_plot} call. The invisible object is a list of matrices, one matrix for each covariate-level/time-slice.
}
\references{
Nason, G.P., Salnikov, D. and Cortina-Borja, M. (2023)
New tools for network time series with an application to COVID-19
hospitalisations. \url{https://arxiv.org/abs/2312.00530}
}

\author{Guy Nason and Daniel Salnikov
%%  ~~who you are~~
}

\examples{

\dontrun{
#
# Produces a R-Corbit plot, which compares three stationary GNAR simulations, where 
# the underlying network is fiveNet.
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
                
# Produce NACF R-Corbit plot with the same network and weights matrix
r_corbit_plot(list(sim1, sim2, sim3), list(fiveNet), 10, 3, list(W), 
		c("sim1", "sim2", "sim3"), same_net = "yes")
# 
# Produce PNACF R-Corbit with different networks and weight matrices
print(r_corbit_plot(list(sim1, sim2, sim3), list(fiveNet, fiveNet, fiveNet), 10, 3, list(W, W, W),
		 c("sim1", "sim2", "sim3"), same_net = "no", partial = "yes"))

}

}


\keyword{ nacf }
\keyword{ pnacf }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ GNAR }
\concept{ R-Corbit }
\concept{ Corbit }
% Use only one concept per line.
