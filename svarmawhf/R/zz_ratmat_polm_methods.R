#' Polynomial Degree
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Compute the polynomial degrees (of the elements) of a polynomial matrix. 
#' Note that for a (scalar) polynomial with zero coefficients the degree 
#' is set to \eqn{-1}.
#' 
#'
#' @param x A polynomial matrix, i.e. an object of class \code{\link{polm}}.
#' @param which (character string) decides whether a matrix with the respectives degrees 
#'              of the entries of the matrix, or a vector with the respective maximal 
#'              degrees in each row or column, or simply the maximum degree of all 
#'              elements of the polynomial matrix is computed. 
#'
#' @return The outcome depends on the parameter \code{which}: 
#' \item{elements}{A matrix with the degrees of the respective elements 
#'                of the polynomial matrix.}
#' \item{rows}{A vector with the maximum degrees within each row.}
#' \item{columns}{A vector with the maximum degrees within each column.}
#' \item{matrix}{maximum of the degrees of the elements of the matrix.}
#'
#' @export
#'
#' @examples
#' x = polm(array(c(0,1,1,0,
#'                  0,0,1,0,
#'                  0,0,0,1,
#'                  0,0,0,0), dim = c(2,2,4)))
#' x
#' degree(x)
#' degree(x, 'rows')
#' degree(x, 'columns')
#' degree(x, 'matrix')
degree = function(x, which = c('elements', 'rows', 'columns', 'matrix')) {
  if (!inherits(x, 'polm')) {
    stop('argument "x" must be a "polm" object!')
  }
  which = match.arg(which)
  x = unclass(x)
  # degree of a univariate polynomial (= vector): 
  # length of vector - 1 - number of zero leading coefficients
  deg_scalar = function(x) {
    length(x) - sum(cumprod(rev(x) == 0)) - 1
  }
  deg = apply(x, MARGIN = c(1,2), FUN = deg_scalar)
  if (which == 'matrix') return(max(deg))
  if (which == 'columns') return(apply(deg, MARGIN = 2, FUN = max))
  if (which == 'rows') return(apply(deg, MARGIN = 1, FUN = max))
  return(deg)
}

#' Column End Matrix of a Polynomial Matrix
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' The \emph{column end matrix} of an \eqn{(m,n)}-dimensional polynomial matrix 
#' \eqn{a(z)=a_0 + a_1 z + \cdots + a_p z^p}{a(z)=a[0] + a[1] z + \dots + a[p] z^p} is defined as follows. 
#' Suppose that the maximum degree of the elements in the \eqn{i}-th column is \eqn{p_i}{p[i]}. Then 
#' the column end matrix is the \eqn{(m,n)} matrix with \eqn{i}-th column equal to the 
#' \eqn{i}-th column of the coefficient matrix \eqn{a_{p_i}}{a[p[i]]}. If a column of 
#' \eqn{a(z)} is zero, then
#' the elements of the corresponding column of the column end matrix are set to \code{NA}'s.
#'
#' @param x A polynomial matrix, i.e. an object of class \code{\link{polm}}.
#'
#' @return The column end matrix.
#' @export
#'
#' @examples
#' x = polm(array(c(0,1,1,0,
#'                  0,0,1,0,
#'                  0,0,0,1,
#'                  0,0,0,0), dim = c(2,2,4)))
#' x
#' degree(x)
#' degree(x, 'columns')
#' col_end_matrix(x)
col_end_matrix = function(x) {
  if (!inherits(x, 'polm')) {
    stop('argument "x" must be a "polm" object!')
  }
  d = dim(x)
  x = unclass(x)
  NAvalue = ifelse(is.complex(x), NA_complex_, NA_real_)
  m = matrix(NAvalue, nrow = d[1], d[2])
  if (length(x) == 0) {
    return(m)
  }

  # degree of a univariate polynomial (= vector): 
  # length of vector - 1 - number of zero leading coefficients
  deg_scalar = function(x) {
    length(x) - sum(cumprod(rev(x) == 0)) - 1
  }
  deg = apply(x, MARGIN = c(1,2), FUN = deg_scalar)
  col_deg = apply(deg, MARGIN = 2, FUN = max)
  
  for (i in iseq(1, dim(x)[2])) {
    if (col_deg[i] >= 0) m[,i] = x[,i,col_deg[i]+1]
  }
  return(m)
}

#' Prune Matrix Polynomial
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Performs three steps to simplify a matrix polynomial.
#' \enumerate{
#' \item All leading coefficients where the absolute values of the real and the imaginary parts are
#'   less than or equal to \code{tol} are set to zero.
#' \item The zero leading coefficient matrices are dropped.
#' \item If all the absolute values of the imaginary parts of the coefficients are less than or
#'       equal to \code{tol} then the coefficients are set to real values.
#' }
#' Empty polynomial matrices (i.e. matrices with zero rows or columns) are set to polynomial of zero degree.
#'
#' @param x \code{\link{polm}} object.
#' @param tol Double. Tolerance paramater. Default set to \code{sqrt(.Machine$double.eps)}.
#' @param brutal Boolean. Default set to FALSE.
#'   If TRUE, all small elements are set to zero (irrespective of whether they are leading
#'   coefficients or not).
#'
#' @return A matrix polynomial, i.e. a \code{\link{polm}} object.
#' @export
#'
#' @examples
#' x = polm(array(c(1,0,0,0,
#'                  0,1,0,0,
#'                  0,0,1,0,
#'                  0,0,0,1,
#'                  0,0,0,0), dim = c(2,2,5)) + 1e-4)
#' x
#' prune(x, tol = 1e-3)
#' prune(x, tol = 1e-3, brutal = TRUE)
#'
#' # Case of complex variables:
#' x = x + complex(imaginary = 1e-5)
#' x
#' prune(x, tol = 1e-3)
#'
#' # also works for constant matrix polynomials (i.e. matrices)
#' x = polm(array(0:3, dim = c(2,2,1))+1e-4)
#' x
#' prune(x, tol = 1e-3)
#' 
#' # empty polynomials are coerced to polynomials of degree zero 
#' x = polm(array(0, dim = c(0,2,5)))
#' x
#' prune(x)
prune = function(x, tol = sqrt(.Machine$double.eps), brutal = FALSE) {
  if (!inherits(x, 'polm')) {
    stop('argument "x" must be an polm object!')
  }
  x = unclass(x)
  d = dim(x)
  if (min(d) <= 0) {
    # empty polynomial, or polynomial of degree (-1)
    return(polm(array(0, dim = c(d[1], d[2], 0))))
  }

  # step one: Set all small leading coefficients to zero
  issmall = ( (abs(Re(x)) <= tol) & (abs(Im(x)) <= tol) )
  issmall = apply(issmall, MARGIN = c(1,2), FUN = function(x) { rev(cumprod(rev(x))) }) 
  # apply: returns an array of dim (d[3], d[1], d[2]) if d[3] > 1
  
  # make sure that issmall is an array (also in the case where the matrix polynomial is constant)
  dim(issmall) = d[c(3,1,2)]
  
  # permute the dimensions back to the polm form: 
  # necessary because apply returns an array of dim (d[3], d[1], d[2]) if d[3] > 1
  issmall = aperm(issmall, c(2,3,1))
  issmall = (issmall == 1)
  
  # finish step one
  x[issmall] = 0
  
  # step two: drop leading zero matrix coefficients
  keep = apply(!issmall, MARGIN = 3, FUN = any)
  # keep[1] = TRUE # keep the constant
  x = x[,, keep, drop = FALSE]
  
  # step three: drop imaginary part if all imaginary parts are small
  if (is.complex(x)) {
    if ( all(abs(Im(x)) <= tol) ) {
      x = Re(x)
    }
  }
  
  # This option is provided to see, e.g., the lower triangularity of the zero power coefficient 
  # matrix when using "transform_lower_triangular"
  if (brutal){
    issmall_brutal = ( (abs(Re(x)) <= tol) & (abs(Im(x)) <= tol) )
    x[issmall_brutal] = 0
  }
  
  x = polm(x)
  return(x)
}


#' Companion Matrix of a Polynomial Matrix
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Computes a companion matrix for a square (\eqn{m,m)}-dimensional), matrix polynomial 
#' \deqn{a(z) = a_0 + a_1 z + \cdots + a_p z^p}{a(z) = a[0] + a[1] z + \dots + a[p] z^p}
#' The companion materix is e.g. used to determine the zeroes of a polynomial matrix,
#' see \code{\link{zeroes}}. 
#' \cr
#' Note that the function throws an error, if the constant term \eqn{a_0}{a[0]} is not regular.
#' There is no check whether some of the leading coefficients are zero. 
#' So the results is an \eqn{(mp,mp)}-dimensional matrix, even if \eqn{a_p}{a[p]} is zero. 
#'
#' @param a A square polynomial matrix, i.e. an object of class \code{\link{polm}}. 
#'
#' @return A (companion) matrix of dimensions \eqn{(mp,mp)}. 
#' @export
#'
#' @examples
#' companion_matrix(polm(c(1,0,0,0.5,0))) # scalar polynomial
#' companion_matrix(polm(diag(3)))        # zero degree polynomial 
#' companion_matrix(polm(dbind(d = 3, diag(2), -test_array(dim = c(2,2,1)))))
#' companion_matrix(polm(dbind(d = 3, diag(2), -test_array(dim = c(2,2,2)))))
#' 
#' \dontrun{
#' # the following examples throw an error
#' companion_matrix(polm(c(0,0,0,0.5))) # constant term is zero
#' companion_matrix(polm(test_array(dim = c(2,1,3)))) # non-square polynomial
#' companion_matrix(polm(test_array(dim = c(2,2,0)))) # zero polynomial
#' }
companion_matrix = function(a) {
  if (!inherits(a, 'polm')) {
    a = try(polm(a), silent = TRUE)
    if (inherits(a, 'try-error')) stop('argument "a" is not coercible to polm object!')
  }
  a = unclass(a)
  d = dim(a)
  if ((d[1] != d[2]) || (d[3] <= 0)) stop('argument "a" must represent a square, non-singular polynomial matrix')
  
  m = d[1]
  p = d[3] - 1
  
  if (m > 0) {
    # check a(0)
    a0 = try(solve(matrix(a[,,1], nrow = m, ncol = m)), silent = TRUE)
    if (inherits(a0, 'try-error')) stop('constant term a[0] is non invertible')
  }

  if ((m*p) == 0) return(matrix(0, nrow = 0, ncol = 0))

  # coerce to (m,m(p+1)) matrix
  dim(a) = c(m, m*(p+1))

  # normalize constant term a[0] -> I, a[i] -> - a[0]^{-1} a[i]
  a = (-a0) %*% a[, (m+1):(m*(p+1)), drop = FALSE]

  if (p == 1) {
    return(a)
  }
  return( rbind(a, diag(x = 1, nrow = m*(p-1), ncol = m*p)) )
}


# internal function 
# univariate polynomial division c = a / b 
poly_div = function(a, b) {
  # a,b are vectors
  
  # take care of the case that the leading coefficients of b are zero ( e.g. b = c(1,2,0,0))
  if (length(b) > 0) {
    b = b[rev(cumprod(rev(b == 0))) == 0]
  }
  lb = length(b)
  if (lb == 0) {
    stop('illegal polynomial division (b is zero)')
  }

  # take care of the case that the leading coefficients of a are zero ( e.g. a = c(1,2,0,0))
  if (length(a) > 0) {
    a = a[rev(cumprod(rev(a == 0))) == 0]
  }
  la = length(a)
  
  if (la < lb) return(0)   # deg(a) < deg(b)
  
  if (lb == 1) return(a/b) # deg(b) = 0
  
  a = rev(a)
  b = rev(b)
  c = rep.int(0, la - lb + 1)
  i = la - lb + 1
  while (i > 0) {
    d = a[1]/b[1]
    c[i] = d
    a[1:lb] = a[1:lb] - d*b
    a = a[-1]
    i = i - 1
  }
  return(c)
}


# internal function 
# univariate polynomial remainder r: a = b * c + r
poly_rem = function(a, b) {
  # a,b are vectors
  
  # take care of the case that the leading coefficients of b are zero ( e.g. b = c(1,2,0,0))
  if (length(b) > 0) {
    b = b[rev(cumprod(rev(b == 0))) == 0]
    lb = length(b)
  } else {
    lb = 0
  }
  if (lb == 0) {
    stop('illegal polynomial division (b is zero)')
  }
  
  # take care of the case that the leading coefficients of a are zero ( e.g. a = c(1,2,0,0))
  if (length(a) > 0) {
    a = a[rev(cumprod(rev(a == 0))) == 0]
    la = length(a)
  } else {
    la = 0
  }

  if (la < lb) return(a)   # deg(a) < deg(b)

  if (lb == 1) {
    return(0)
  }
  
  a = rev(a)
  b = rev(b)
  while (length(a) >= lb) {
    d = a[1]/b[1]
    a[1:lb] = a[1:lb] - d*b
    a = a[-1]
  }
  return( rev(a) )
}

# internal function
# l2 norm of a vector
l2_norm = function(x) return(sqrt(sum(x^2)))

# internal function
# consider a vector x = c(x[1], ..., x[k],  x[k+1], ..., x[n])
# return a logical  i = c(FALSE, .., FALSE, TRUE, .., TRuE) 
# where k is the minimum integer such that | x[s] | <= tol for all s > k 
is_small = function(x, tol = sqrt(.Machine$double.eps), count = TRUE) {
  if (length(x) == 0) {
    i = logical(0)
  } else {
    i = rev( cumprod( rev(abs(x) <= tol) ) == 1 )
  }
  if (count) {
    return(sum(i))
  } else {
    return(i)
  }
}

# internal function
# consider a vector x = c(x[1], ..., x[k],  x[k+1], ..., x[n])
# return a logical  i = c(TRUE, .., TRUE, FALSE, .., FALSE) 
# where k is the minimum integer such that | x[s] | <= tol for all s > k 
is_large = function(x, tol = sqrt(.Machine$double.eps), count = TRUE) {
  if (length(x) == 0) {
    i = logical(0)
  } else {
    i = rev( cumprod( rev(abs(x) <= tol) ) == 0 )
  }
  if (count) {
    return(sum(i))
  } else {
    return(i)
  }
}