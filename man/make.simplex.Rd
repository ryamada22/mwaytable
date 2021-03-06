% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SphereKaji20140527.R
\name{make.simplex}
\alias{make.simplex}
\title{Make a rotation matrix using simple vertices coordinates (2)}
\usage{
make.simplex(k)
}
\arguments{
\item{k}{A integer}
}
\value{
k x k rotation matrix to transfer k vertices so that all k vertices' k-th coordinates are the same
}
\description{
Make a rotation matrix using simple vertices coordinates (2)
}
\examples{
k <- 3
M <- make.simplex(k)
M \%*\% diag(rep(1,k))
}

