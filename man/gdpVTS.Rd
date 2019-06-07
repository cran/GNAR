\name{gdpVTS}
\alias{gdpVTS}
\title{Differenced GDP values for 35 countries}
\usage{
gdpVTS
}
\description{
This dataset is from the OECD (OECD (2018), Quarterly GDP (indicator). <doi:10.1787/b86d1fc8-en> (Accessed on 29 January 2018)) and is differenced annual growth rate for 35 countries for 1962-2013.
}
\format{
\code{gdpVTS} is a \link[stats]{ts} object with a matrix of 52 rows (t=52) and 35 columns (n=35)
}
\examples{
#Plot using 'ts' S3 function, can only plot up to 10 columns at once
plot(gdpVTS[,1:5])

#Plot as heatmap
image(gdpVTS)
}
