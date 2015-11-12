% Generated by roxygen2 (4.0.1): do not edit by hand
\name{map2df}
\alias{map2df}
\title{Make a matrix that transfer a vectorized table array to a vector in the df-dimensional space}
\usage{
map2df(A)
}
\arguments{
\item{A}{An array or matrix of multiway table}
}
\value{
A matrix
}
\description{
Make a matrix that transfer a vectorized table array to a vector in the df-dimensional space
}
\examples{
r <- c(2,3,4)
A <- array(1:prod(r),r)
M <- map2df(A)
M \%*\% c(A)
}
