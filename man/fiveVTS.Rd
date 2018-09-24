\name{fiveVTS}
\alias{fiveVTS}
\title{Simulated GNAR data using fiveNet}
\usage{
fiveVTS
}
\description{
This dataset is simulated from the GNAR(1,[1]) model, using fiveNet as the accompanying network. The parameters used were \eqn{\alpha_1 = 0.2} for each node, and \eqn{\beta_1} = 0.5. See \link[GNAR]{GNARfit} for the model description.
}
\format{
\code{fiveVTS} is a \link[stats]{ts} object with a matrix of 200 rows (t=200) and 5 columns (n=5)
}
\examples{
data(fiveNode)
#Plot using ts S3 function
plot(fiveVTS)

#Plot as heatmap
image(fiveVTS)
}
