\name{logMVbedMVC.vts}
\docType{data}
\alias{logMVbedMVC.vts}
\title{Number of COVID19 admission to mechanical ventilation beds}
\description{
This matrix contains a multivariate/vector time series that counts
the number of daily admissions to mechanical ventilation beds in one of
140 NHS Trusts from 2nd April 2020 to 27th June 2021.
}
\format{
The dimension
of the matrix is 452x140. The dates and Trust ID codes are stored in
the dimnames (first and second) respectively.
}
\source{UK Coronavirus website \url{https://ukhsa-dashboard.data.gov.uk} }
\references{
Nason, G.P., Salnikov, D. and Cortina-Borja, M. (2023)
New tools for network time series with an application to COVID-19
hospitalisations. \url{https://arxiv.org/abs/2312.00530}
}
\section{SEE ALSO}{\code{\link{corbit_plot}},\code{\link{NHSTrustMVCAug120.net}}}
\examples{

\dontrun{
data(logMVbedMVC.vts)
data(NHSTrustMVCAug120.net)

#
# Do a corbit plot with this data, with only three lags and one stage
#
# Note, normally max_lag and max_stage would be higher than this, the
# values are artificially small here, as otherwise the run-time restrictions
# for CRAN packaging might be exceeded. 

corbit_plot(vts=logMVbedMVC.vts, net=NHSTrustMVCAug120.net, max_lag=3, max_stage=1)
}

}
\keyword{datasets}
\author{Daniel Salnikov and Guy Nason}
