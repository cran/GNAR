\name{vswind}
\alias{vswind}
\alias{vswindcoords}
\alias{vswindnames}
\alias{vswindnet}
\alias{vswindnetD}
\alias{vswindts}

\docType{data}
\title{
Wind Speed example network time series
}
\description{A suite of data objects concerning wind speed analysis.
	The dataset contains a multivariate time series of wind speeds,
	two network descriptions, a vector of names for weather stations,
	and the coordinates of the weather stations.
}
\usage{data("vswind")}
\format{ This dataset contains six R objects:\cr
\code{vswindts} is a \link[stats]{ts} object with a matrix of 721 rows (t=721) and 102 columns (n=102). This corresponds to 721 observations made through time at
102 weather stations.
\code{vswindnetD} is a \code{GNARnet}
	object containing \code{$edges} and \code{$dist}. \cr
\code{edges} is a list of length 102, with \code{edges[[i]]} containing the vertices that node i is connected to. \cr
\code{dist} is a list of length 102, with \code{dist[[i]]} containing the length of the vertices that node i is connected to.
\code{vswindnet} is the same as \code{vswindnetD} except that all the distances
are replaced by 1.
\code{vswindnames} is a character vector of length 102 containing the wind speed
	site names and
\code{vswindcoords} is a matrix with 102 rows (one for each wind station) and
two columns providing the x and y coordinates of the weather stations.
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
	The base data were obtained from the
	\url{http://wow.metoffice.gov.uk} UK Met Office WeatherObservationsWebsite
	distributed under the UK Open Government License
	\url{https://www.nationalarchives.gov.uk/doc/open-government-licence/version/1/open-government-licence.htm}
	Contains public sector information licensed under the Open Goverment
	Licence v1.0.

}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\seealso{\code{\link{windnetplot}}}
\examples{
#
# The name entry for Bristol
#
vswindnames[77]
#[1] "BRIST"
#
# plot the distance network
#
\dontrun{windnetplot()}
}
\keyword{datasets}
