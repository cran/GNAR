\name{na.row}
\alias{na.row}
\title{Identifies which rows of a matrix have NAs}
\usage{
na.row(mat)
}
\description{
Returns a vector with elements TRUE/FALSE identifying which rows contain NA elements.
}
\arguments{
\item{mat}{a matrix object.}
}
\details{This function is used in the unstacking of residuals into a residual matrix and replacing NAs where they were previously present.}
\value{
\code{na.row} returns a vector of length equal to the number of rows in \code{mat}. Each element is either TRUE or FALSE.
}
\examples{
#Check if there are and NAs in fiveVTS
na.row(fiveVTS)
}
