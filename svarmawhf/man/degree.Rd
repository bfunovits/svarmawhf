% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zz_ratmat_polm_methods.R
\name{degree}
\alias{degree}
\title{Polynomial Degree}
\usage{
degree(x, which = c("elements", "rows", "columns", "matrix"))
}
\arguments{
\item{x}{A polynomial matrix, i.e. an object of class \code{\link{polm}}.}

\item{which}{(character string) decides whether a matrix with the respectives degrees 
of the entries of the matrix, or a vector with the respective maximal 
degrees in each row or column, or simply the maximum degree of all 
elements of the polynomial matrix is computed.}
}
\value{
The outcome depends on the parameter \code{which}: 
\item{elements}{A matrix with the degrees of the respective elements 
               of the polynomial matrix.}
\item{rows}{A vector with the maximum degrees within each row.}
\item{columns}{A vector with the maximum degrees within each column.}
\item{matrix}{maximum of the degrees of the elements of the matrix.}
}
\description{
This function was originally part of the R-package \strong{rationalmatrices}.
\cr
Compute the polynomial degrees (of the elements) of a polynomial matrix. 
Note that for a (scalar) polynomial with zero coefficients the degree 
is set to \eqn{-1}.
}
\examples{
x = polm(array(c(0,1,1,0,
                 0,0,1,0,
                 0,0,0,1,
                 0,0,0,0), dim = c(2,2,4)))
x
degree(x)
degree(x, 'rows')
degree(x, 'columns')
degree(x, 'matrix')
}
