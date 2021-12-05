# poles.___ and zeroes.___ method ############################################################
# 
# 
#' Poles and Zeroes
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Compute the poles and zeroes of a rational matrix. For polynomial matrices and rational matrices 
#' in left matrix fraction form the poles (and zeroes) are computed via the (reciprocals of the) 
#' eigenvalues of the associated companion matrices, see also \code{\link{companion_matrix}}. 
#' For statespace realizations the poles are computed via (the reciprocals of) the 
#' eigenvalues of the state transition matrix \eqn{A} and for the zeroes 
#' the eigenvalues of the state transition matrix \eqn{A-BD^{-1}C} of the inverse are used. 
#' 
#' The methods do not return numerically reliable and correct results in all cases. 
#' For some more details see the vignette \href{../doc/rational_matrices.html}{Rational Matrices}. 
#' 
#' \itemize{
#' \item Zeroes are only computed for square, non singular matrices which have no zero 
#'       at \eqn{z=0}. If the matrix evaluated at \eqn{z=0} is close to singular, 
#'       the results may be unreliable. 
#' \item The procedures use a threshold \code{tol} in order to decide whether 
#'       a small eigenvalue returned by \code{\link{eigen}} corresponds to 
#'       a "true zero" eigenvalue or not. 
#' \item If the pair \eqn{a,b} of polynomials of the LMFD is not left coprime then 
#'       a pole/zero cancellation occurs. This is not taken into account by the procedures. Hence, 
#'       in this case, the results also contain some spurious poles/zeroes. This happens 
#'       also for non minimal state space realizations.      
#' }
#'
#' @param x an object which represents a rational matrix 
#'          (i.e. a \code{\link{polm}}, \code{\link{lmfd}} or \code{\link{stsp}} object).
#' @param tol Double. Default set to \code{sqrt(.Machine$double.eps)}.
#'   Required to decide on when a root is to be considered "at infinity".
#' @param print_message Boolean. Default set to TRUE.
#'   Prints a message if roots "at infinity " are discarded.
#' @param ... not used.
#'
#' @return Vector of poles, respectively zeroes.
#'
#' @rdname poles_and_zeroes
#' @name poles and zeroes
#' @export
poles = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  UseMethod("poles", x)
}

#' @rdname poles_and_zeroes
#' @export
zeroes = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  UseMethod("zeroes", x)
}


#' @rdname poles_and_zeroes
#' @export
#' 
#' @examples
#'  
#' # zeroes of polynomial matrices #############################################
#' 
#' # scalar polynomial ###
#' a = polm(c(1, 0, 0, 0.5, 0))
#' (z = zeroes(a))
#' 
#' # compare with the result of "polyroot"
#' all.equal(sort(z), sort(polyroot(as.vector(a))))
#' 
#' # zero degree polynomial (have no zeroes) ###
#' zeroes(polm(diag(3)))
#' 
#' # (2 x 2) polynomial of degree 2 ### 
#' a = polm(dbind(d = 3, diag(2), test_array(dim = c(2,2,2))))
#' (z = zeroes(a))
#' 
#' \dontrun{
#' # the following examples throw an error
#' zeroes(polm(c(0, 0, 0, 0.5))) # constant term is zero
#' zeroes(polm(test_array(dim = c(2, 1, 3)))) # non-square polynomial
#' zeroes(polm(test_array(dim = c(2, 2, 0)))) # zero polynomial
#' }
zeroes.polm = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  x = prune(x, tol = 0) # skip zero coefficients
  d = dim(unclass(x))
  if ((d[1] != d[2]) || (d[1] == 0) || (d[3] <= 0)) {
    stop('argument "x" must represent a non-empty, square and non-singular polynomial matrix')
  }
  
  # compute companion matrix
  A = try(companion_matrix(x), silent = TRUE)
  if (inherits(A,'try-error')) {
    stop('could not generate companion matrix (constant term a[0] is singular?)')
  }
  
  # zero degree polynomial
  if (nrow(A) == 0) return(numeric(0))

  z = eigen(A, only.values=TRUE)$values
  
  if (any(abs(z) <= tol)){
    z <- z[!(abs(z) <= tol)]
    if (print_message){
      message("There are determinantal roots at (or close to) infinity.\nRoots close to infinity got discarded.")
    }
  }
  
  return(1/z)
}


#' @rdname poles_and_zeroes
#' @export
#' 
#' @examples
#' 
#' # poles of polynomial matrices #############################################
#' 
#' # polynomials have no poles ###
#' poles(test_polm(dim = c(2,1), degree = 2, random = TRUE)) 
poles.polm = function(x, ...) {
  # polynomial matrices have no poles 
  return(numeric(0))
}


#' @rdname poles_and_zeroes
#' @export
#' 
#' @examples
#'  
#' # zeroes of a rational matrix in LMFD form #################################
#' 
#' c = lmfd(test_polm(dim = c(2,2), degree = 3, random = TRUE),
#'          test_polm(dim = c(2,2), degree = 1, random = TRUE))
#' (z = zeroes(c))
#' all.equal(z, zeroes(c$b))
zeroes.lmfd = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  d = dim(x)
  if ((d[1] != d[2]) || (d[1] == 0)) {
    stop('argument "x" must represent a non-empty, square and non-singular rational matrix in LMFD form')
  }
  return(zeroes(x$b, tol, print_message))
}

#' @rdname poles_and_zeroes
#' @export
#' 
#' @examples
#' 
#' # poles of a rational matrix in LMFD form ##################################
#' 
#' (z = poles(c))
#' all.equal(z, zeroes(c$a))
poles.lmfd = function(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...) {
  d = dim(x)
  if ((d[1] == 0) || (d[3] < 0)) {
    stop('argument "x" does not represent a valid LMFD form')
  }
  return(zeroes(x$a, tol, print_message))
}
