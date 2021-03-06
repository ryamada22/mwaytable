% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SphereKaji20140527.R
\name{make.test.vecs}
\alias{make.test.vecs}
\title{Generate test vectors from tables}
\usage{
make.test.vecs(t.sphere, tests)
}
\arguments{
\item{t.sphere}{table.sphere function's output}

\item{tests}{A list of arrays/matrix of tables representing tests}
}
\value{
A number-of-test rows x df columns matrix of test vectos, each row of which is a test vector

alt.vec The vector of the estimated table in df-dimensional space
}
\description{
Generate test vectors from tables
}
\examples{
A <- matrix(c(10,20,30,40,50,60),2,3)
t.sphere <- table.sphere(A)
t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
t2[1,2] <- 0.5
t3[1,2] <- 1
tests <- list(t1,t2,t3)
make.test.vecs(t.sphere,tests)
}

