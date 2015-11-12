% Generated by roxygen2 (4.0.1): do not edit by hand
\name{perm.data}
\alias{perm.data}
\title{Randomize a data matrix}
\usage{
perm.data(data, v)
}
\arguments{
\item{data}{a matrix of sample size x number of variables}

\item{v}{a vector of 0 or 1 with length being the number of variables, which indicates variables to be shuffled.}
}
\value{
data matrix
}
\description{
Randomize a data matrix
}
\examples{
m <- matrix(c(10,20,30,40,50,60),2,3)
d <- make.data.from.table(m)
make.table.from.data(d)
}
