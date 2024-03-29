\name{GNAR}
\alias{GNAR}

\title{GNAR package}
\description{A package to fit, predict, and simulate time series using the Generalised Network AutoRegressive (GNAR) model, potentially with exogenous variables. The main functions are \link[GNAR]{GNARfit} and \link[GNAR]{GNARXfit}, which fits the model to a time series and network(s), S3 method \link[GNAR]{predict.GNARfit} which predicts from a fitted GNAR model, and \link[GNAR]{GNARsim} which simulates from a GNAR model with specified parameters. For details of the model, see \link[GNAR]{GNARfit}. The package also contains an example network time series in data file \link[GNAR]{fiveNode}, with network \link[GNAR]{fiveNet}, and simulated time series \link[GNAR]{fiveVTS}.}
\references{
Knight, M.I., Nunes, M.A. and Nason, G.P. (2015) Modelling, detrending and decorrelation of network time series. \href{https://arxiv.org/abs/1603.03221}{arXiv preprint}.\cr\cr
Knight, M.I., Leeming, K., Nason, G.P. and Nunes, M. A. (2020) Generalised Network Autoregressive Processes and the GNAR package. \emph{Journal of Statistical Software}, \bold{96} (5), 1--36.\cr\cr
Nason G.P. and Wei J. (2022) Quantifying the economic response to COVID-19 mitigations and death rates via forecasting Purchasing Managers’ Indices using Generalised Network Autoregressive models with exogenous variables. \emph{Journal of the Royal Statistical Society Series A}, \bold{185} (4), 1778--1792.
}
