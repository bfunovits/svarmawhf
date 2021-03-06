% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zz_ratmat_dim_methods.R
\name{dim methods}
\alias{dim methods}
\alias{dim.polm}
\alias{dim.lmfd}
\alias{dim.pseries}
\title{Dimensions of Objects}
\usage{
\method{dim}{polm}(x)

\method{dim}{lmfd}(x)

\method{dim}{pseries}(x)
}
\arguments{
\item{x}{object}
}
\value{
returns a named integer vector
\item{\code{c(m,n,p)}}{for polynomial matrices (i.e. for \code{\link{polm}} objects).}
\item{\code{c(m,n,p,q)}}{for left matrix fraction descriptions (i.e. for \code{\link{lmfd}} objects).}
\item{\code{c(m,n,s)}}{for state space representations (i.e. for \code{\link{stsp}} objects).}
\item{\code{c(m,n,lag.max)}}{for impulse response functions (i.e. for \code{\link{pseries}} objects).}
\item{\code{c(m,n,n.f)}}{for frequency response functions (i.e. for \code{\link{zvalues}} objects).}
}
\description{
This function was originally part of the R-package \strong{rationalmatrices}.
}
