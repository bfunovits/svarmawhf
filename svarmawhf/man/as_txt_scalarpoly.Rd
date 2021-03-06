% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zz_rldm_print_methods.R
\name{as_txt_scalarpoly}
\alias{as_txt_scalarpoly}
\title{Coerce Scalar Polynomials to Character Strings}
\usage{
as_txt_scalarpoly(coefs, syntax = c("txt", "TeX", "expression"), x = "z")
}
\arguments{
\item{coefs}{(numeric) vector of coefficients.}

\item{syntax}{(character string) determines the format of the output string.}

\item{x}{(character string) polynomial variable.}
}
\value{
character string.
}
\description{
This function was originally part of the R-package \strong{RLDM}.
\cr
This utility coerces a scalar polynomial (given by the vector of coefficients)
to a character string. The following "formats" are implemented.
\code{syntax = "txt"} returns a simple text representation,
\code{syntax = "TeX"} renders the coefficients to string in "TeX" syntax and
\code{syntax = "expression"} gives a string which may be rendered to
an \code{R} expression with \code{\link[base]{parse}}. This expression
may be used to evaluate the polynomial and for annotating plots,
see \code{\link[grDevices]{plotmath}} and the examples below.
}
\examples{
coefs = c(1, 2.3, 0, -1, 0)

as_txt_scalarpoly(coefs, syntax = 'txt', x = 'x')
as_txt_scalarpoly(coefs, syntax = 'TeX', x = '\\\\alpha')
as_txt_scalarpoly(coefs, syntax = 'expression', x = 'z')

\dontrun{
# the case syntax = "expression" may be used e.g. as follows

# make_polyfun creates a "closure" which evaluates the polynomial at given points
# note that this simple version does not work for zero polnomials!
make_polyfun = function(coefs) {
  expr = parse(text = as_txt_scalarpoly(coefs, 'expression', 'x'))
  fun = function(x) {
    return(eval(expr))
  }
  return(fun)
}

a = make_polyfun(coefs)
a(1)   # return the value  of the polynomial at x = 1
a(1:5) # return the values of the polynomial at x = 1,2,3,4,5

# create a plot
x_grid = seq(from = -1, to = 1, length.out = 101)
plot(x_grid, a(x_grid), type = 'l', xlab = 'x', ylab = 'a(x)',
     main = parse(text = paste('a(x) == ',
          as_txt_scalarpoly(coefs, syntax = 'expression', x = 'x'))))
}
}
