# classes.R ##########################################################################################

#' Constructor for Polynomial Matrices
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' \code{polm} objects represent polynomial matrices
#' \deqn{a(z) = a_0 + a_1 z + \cdots + a_p z^p}{a(z) = a[0] + a[1]z + \dots + a[p]z^p}
#' If the matrix \eqn{a(z)} is an \eqn{(m,n)} dimensional polynomial matrix (i.e. the coefficients 
#' \eqn{a_i}{a[i]} are \eqn{(m,n)} dimensional real or complex valued matrices) then the 
#' \code{polm} object stores the coefficients in an \code{(m,n,p+1)} dimensional (numeric or complex) 
#' array together with a class attribute \code{c("polm","ratm")}. 
#' \cr
#' The constructor function \code{polm(a)} takes a (numeric or complex) vector, matrix or 
#' 3-dimensional array and returns a \code{polm} object.
#' 
#' Any of the dimensions of the 3-dimensional array may also be zero. In particular, 
#' if the third dimension is zero, then the \code{polm} object is interpreted as the zero polynomial.
#' 
#' For important methods and functions for this class have a look at the "see also" section.
#' 
#' @param a either a (numeric or complex) vector, matrix or 3-D array. A vector is coerced to a scalar 
#'          (i.e. \eqn{(1,1)}-dimensional) polynomial and a matrix gives a polynomial matrix of zero degree. 
#'
#' @return An object of class \code{polm}. 
#' @export
#' 
#' @seealso 
#' \itemize{
#'    \item \code{\link{test_polm}} generates random polynomials. 
#'    \item checks: \code{\link{is.polm}}, \code{\link{is.miniphase}}, 
#'          and \code{\link{is.coprime}}. As a byproduct \code{is.coprime(a)} 
#'          computes the zeroes of a square polynomial matrix \eqn{a(z)}. 
#'    \item generic S3 methods: \code{\link{dim}}, 
#'          \code{\link{str}} and \code{\link{print}}.
#'    \item arithmetics: \code{\link{Ops.ratm}}, matrix multiplication \code{\link{\%r\%}}, 
#'          polynomial division \code{\%/\%}, 
#'          polynomial remainder \code{\%\%}, ...
#'    \item matrix operations: \code{\link{t.polm}}, 
#'          \code{\link{bind}}, \code{\link{[.polm}}, \code{\link{[<-.polm}}, ...
#'    \item \code{\link{degree}} returns the degree, \code{\link{col_end_matrix}} computes the 
#'          \emph{column end matrix} and \code{\link{prune}} "simplifies" a polynomial. 
#'    \item \code{\link{blaschkerize}} and \code{\link{transform_allpass}} may be used 
#'          to reflect zeroes of a polynomial matrix by multiplication
#'          with allpass rational matrices.
#'    \item normal forms: Hermite normal form \code{\link{hnf}}, Smith normal form \code{\link{snf}},  
#'                        column reduced form \code{\link{col_reduce}} and Wiener Hopf factorization 
#'                        \code{\link{whf}}.
#'    \item \code{\link{companion_matrix}}, \code{\link{zeroes}}, \code{\link{pseries}}, \code{\link{zvalues}}, ...  
#' }
#' 
#' @examples 
#' # (1 x 1) polynomial matrix a(z) =  0 + 1z + 2z^2
#' polm(0:2)
#' 
#' # (2 x 3) polynomial matrix a(z) = a0 (degree is zero)
#' polm(diag(1, nrow = 2, ncol = 3))
#' 
#' # random (2 x 3) polynomial matrix a(z) = a0 + a1 z + a2 z^2 + a3 z^3 (degree = 3)
#' polm(array(stats::rnorm(2*3*4), dim = c(2,3,4)))
#' 
#' # random (2 x 1) polynomial matrix with complex coefficients (degree = 2)
#' a = polm(array(complex(real = stats::rnorm(2*1*3), 
#'                imaginary = stats::rnorm(2*1*3)), dim = c(2,1,3)))
#' is.polm(a)
#' dim(a)
#' str(a)
#' print(a, digits = 3)
polm = function(a) {
  
  if (!(is.numeric(a) || is.complex(a))) {
    stop('input "a" must be a numeric or complex vector/matrix/array!')
  }
  if (is.vector(a)) {
    dim(a) = c(1,1,length(a))
  }
  if (is.matrix(a)) {
    dim(a) = c(dim(a),1)
  }
  
  if ( (!is.array(a)) || (length(dim(a)) !=3) ) {
    stop('could not coerce input parameter "a" to a valid 3-D array!')
  }
  class(a) = c('polm','ratm')
  return(a)
}

#' Constructor for Left Matrix Fraction Descriptions (LMFDs)
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' A Left Matrix Fraction Description (LMFD) of a rational matrix, \eqn{x(z)} say, is a
#' pair \eqn{(a(z),b(z))} of polynomial matrices, such
#' that \eqn{x(z) = a^{-1}(z) b(z)}. 
#' The polynomial matrix \eqn{a(z)} must be square and invertible.
#' 
#' Suppose that \eqn{x(z)=a^{-1}(z) b(z)} is an \eqn{(m,n)}-dimensional matrix and that 
#' \eqn{a(z)} and \eqn{b(z)} have degrees \eqn{p} and \eqn{q} respectively. 
#' The corresponding \code{lmfd} object stores the coefficients of the polynomials \eqn{a(z), b(z)} in 
#' an \eqn{(m,m(p+1)+n(q+1))} dimensional (real or complex valued) matrix together with an 
#' attribute \code{order = c(m,n,p,q)} and a class attribute \code{c("lmfd", "ratm")}. 
#' 
#' For a valid LMFD we require \eqn{m>0} and \eqn{p\geq 0}{p \ge 0}. 
#'
#' @param a,b \code{\link{polm}} objects, or objects which may be coerced to a \code{\link{polm}} object,
#'            via \code{x = polm(x)}. Either of the two arguments may be omitted.
#'
#'
#' @seealso Useful methods and functions for the \code{lmfd} class are: 
#' \itemize{
#'    \item \code{\link{test_lmfd}} generates random rational matrices in LMFD form.
#'    \item checks: \code{\link{is.lmfd}}, \code{\link{is.miniphase}}, \code{\link{is.stable}} 
#'          and \code{\link{is.coprime}}.
#'    \item generic S3 methods: \code{\link{dim}}, 
#'          \code{\link{str}} and \code{\link{print}}.
#'    \item arithmetics: \code{\link{Ops.ratm}}.
#'    \item matrix operations: \code{\link{bind}}.
#'    \item extract the factors \eqn{a(z)} and \eqn{b(z)} with \code{\link{$.lmfd}}.      
#'    \item \code{\link{zeroes}}, \code{\link{pseries}}, \code{\link{zvalues}}, ...  
#' }
#' 
#' @return An object of class \code{lmfd}.
#' @export
#' 
#' @examples 
#' ### (1 x 1) rational matrix x(z) = (1+z+z^2)^(-1) (3+2z+z^2)
#' lmfd(c(1,1,1), c(3,2,1)) %>% print(format = 'c')
#' 
#' ### (1 x 1) rational matrix x(z) = (3+2z+z^2)
#' lmfd(b = c(3,2,1)) %>% print(format = 'c')
#' 
#' ### (1 x 1) rational matrix x(z) = (1+z+z^2)^(-1)
#' lmfd(c(1,1,1)) %>% print(format = 'c')
#' 
#' ### (2 x 3) rational matrix with degrees p=1, q=1
#' x = lmfd(array(rnorm(2*2*2), dim = c(2,2,2)), 
#'          array(rnorm(2*3*2), dim = c(2,3,2)))
#' is.lmfd(x)
#' dim(x)
#' str(x)
#' print(x, digits = 2)
#' 
#' \dontrun{
#' ### the following calls to lmfd() throw an error 
#' lmfd() # no arguments!
#' lmfd(a = test_polm(dim = c(2,3), degree = 1))  # a(z) must be square 
#' lmfd(a = test_polm(dim = c(2,2), degree = -1)) # a(z) must have degree >= 0
#' }
lmfd = function(a, b) {
  if (missing(a) && missing(b)) {
    stop('no arguments have been provided')
  }
  if (!missing(a)) {
    if (!inherits(a,'polm')) {
      a = try(polm(a))
      if (inherits(a, 'try-error')) {
        stop('could not coerce "a" to a "polm" object!')
      }
    }
    dim_a = dim(unclass(a))
    if ((dim_a[1] == 0) || (dim_a[1] != dim_a[2]) || (dim_a[3] == 0)) {
      stop('"a" must represent a square polynomial matrix with degree >= 0')
    } 
  }
  if (!missing(b)) {
    if (!inherits(b,'polm')) {
      b = try(polm(b))
      if (inherits(b, 'try-error')) {
        stop('could not coerce "b" to a "polm" object!')
      }
    }
    dim_b = dim(unclass(b))
  }
  if (missing(b)) {
    b = polm(diag(dim_a[2]))
    dim_b = dim(unclass(b))
  }
  if (missing(a)) {
    a = polm(diag(dim_b[1]))
    dim_a = dim(unclass(a))
  }
  if (dim_a[2] != dim_b[1]) {
    stop('the arguments "a", "b" are not compatible')
  }
  
  c = matrix(0, nrow = dim_b[1], ncol = dim_a[2]*dim_a[3] + dim_b[2]*dim_b[3])
  if ((dim_a[2]*dim_a[3]) > 0) c[, 1:(dim_a[2]*dim_a[3])] = a
  if ((dim_b[2]*dim_b[3]) > 0) c[, (dim_a[2]*dim_a[3]+1):(dim_a[2]*dim_a[3] + dim_b[2]*dim_b[3])] = b
  
  c = structure(c, order = as.integer(c(dim_b[1], dim_b[2], dim_a[3]-1, dim_b[3] - 1)),  
                class = c('lmfd','ratm'))
  return(c)
}
