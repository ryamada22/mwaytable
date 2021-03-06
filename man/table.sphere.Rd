% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SphereKaji20140527.R
\name{table.sphere}
\alias{table.sphere}
\title{Spherization of multiway table array}
\usage{
table.sphere(A)
}
\arguments{
\item{A}{An array or matrix of multiway table}
}
\value{
tables A list of three tables: data :Input table A itself, expected : its expected table and diff: table of differecen between A and expected table

r.vec A integer vector indicating numbers of levels of A

df Degrees of freedom of A

matrices A list of four matrices, X,X.sub,X.inv and X.inv.sub, which are output of make.X.vector0, rotating to/from number-of-cell dimension and df dimension with/without meaningless rows/columns

zero.vec A vector indicating meaningful rows of rotations

map.matrices A list of four matrices: Full2df is a matrix transforming a table array vector to df-dimensional point, DF2full is a matrix transforming df-dimensional point to a table vector, Normalvec2DF is a matrix transforming test table to df-dimensional vector and W is a key matrix to calculate other matrices
}
\description{
Spherization of multiway table array
}
\examples{
r <- c(2,3,4)
A <- array(1:prod(r),r)
t.sphere <- table.sphere(A)
}

