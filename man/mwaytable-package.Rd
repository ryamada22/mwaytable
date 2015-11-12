\name{mwaytable-package}
\alias{mwaytable-package}
\alias{mwaytable}
\docType{package}
\title{
What the package does (short line)
Probability density of maximum statistic among multiple df=1 statistics on multi-way tables
}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab mwaytable\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-27\cr
License: \tab What license is it under?\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
~~ directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
A <- matrix(c(10,20,30,40,50,60),2,3)
t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
t2[1,2] <- 0.5
t3[1,2] <- 1
tests <- list(t1,t2,t3)
n <- 1000
mway.test(A,tests,n=10000)
}
