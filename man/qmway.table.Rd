% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SphereKaji20140527.R
\name{qmway.table}
\alias{qmway.table}
\title{Estimate quantiles of maximum statistics on multi-way table under unll or alternative hypothesis}
\usage{
qmway.table(p, A, tests, lower.tail, x, nc, one.side, n)
}
\arguments{
\item{p}{probability vector}

\item{A}{An array or matrix of table that specifies marginal counts}

\item{tests}{A list of tables}

\item{lower.tail}{logical; if TRUE (default), probabilities are P[X ?????? x], otherwise, P[X > x].}

\item{x}{a vector of K statistic values}

\item{nc}{Optional. if FALSE (default), random tables are in null hypothesis and otherwise in alternative hypothesis with its ccenter being table A}

\item{one.side}{Optional. When TRUE, tests are one-sided, otherwise two-sided}

\item{n}{number of random unit vectors used for estimation}
}
\value{
\itemize{
\item{"x"}{quantile}
\item{"p.used"}{Argument p}
\item{"p.estimated"}{p value corresponding to the estimated quantile, which is some deviated from p of argument because, estimation process is involeved.}
}
}
\description{
Estimate quantiles of maximum statistics on multi-way table under unll or alternative hypothesis
}
\examples{
A <- matrix(c(10,20,30,40,50,60),2,3)
t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
t2[1,2] <- 0.5
t3[1,2] <- 1
tests <- list(t1,t2,t3)
n <- 1000
p <- 0.05
qmway.table(p,A,tests)
}

