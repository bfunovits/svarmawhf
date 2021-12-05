# tools.R ###################################################
#

#' Create Test ARMA model
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' \cr
#' This simple tool may be used to create a random ARMA model
#' \deqn{y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p} = b_0 u_t + b_1 u_{t-1} + \cdots + b_q u_{t-q}}{
#'       y[t] + a[1] y[t-1] + ... + a[p] y[t-p] = b[0] u[t] + b[1] u[t-1] + ... + b[q] u[t-q]}
#' with given order \eqn{(p,q)}.
#'
#' We require \eqn{m>0} and \eqn{p\geq 0}{p\ge 0}.
#'
#' The \eqn{b_0}{b[0]} matrix defaults to a \eqn{(m,n)}-dimensional diagonal matrix
#' with ones on the diagonal (\code{diag(x=1, nrow = m, ncol = n)}). However, one may
#' also pass an arbitray (compatible) matrix to the procedure.
#' This matrix may contain \code{NA}'s, which then are replaced by random numbers.
#'
#' The \eqn{sigma_L} matrix defaults to a \eqn{(n,n)}-dimensional lower, triangular matrix
#' However, one may also pass an arbitray (compatible) \eqn{sigma_L} matrix to the procedure.
#'
#' The user may prescribe lower bounds for the moduli of the zeroes and/or poles of the transfer function
#' \deqn{k(z) = a^{-1}(z) b(z).}
#' In this case the procedure simply generates (up to n.trials) random models until a model is found
#' which satisfies the constraint. The standard deviation of the normal distribution, which is used to
#' generate the random entries, is decreased in each step. Of course this is a very crude method and
#' it may fail or need a very large number of randomly generated matrices.
#'
#' @param dim integer vector \code{c(m,n)}.
#' @param degrees integer vector \code{c(p,q)}.
#' @param digits integer, if non NULL then the randomly generated numbers are rounded to
#'               "digits" number of decimal places.
#' @param b0     \eqn{(m,n)} dimensional matrix (or \code{NULL}). See the details below.
#' @param sigma_L     \eqn{(n,n)} dimensional matrix (or \code{NULL}). See the details below.
#' @param bpoles lower bound for the moduli of the poles of the corresponding transfer function (or NULL).
#' @param bzeroes lower bound for the moduli of the zeroes of the corresponding tranmsfer function (or NULL).
#'                This parameter is ignored for non-square matrices (m != n).
#' @param n.trials maximum number of trials.
#'
#' @return \code{\link{armamod}} object, which represents the generated ARMA model.
#' @export
#'
#' @examples
#' ### generate a random ARMA(1,1) model (with two outputs)
#' ### we require that the model is stable and minimum phase
#' model = try(test_armamod(dim = c(2,2), degrees = c(1,1), digits = 2, bpoles = 1, bzeroes = 1))
#' if (!inherits(model, 'try-error')) {
#'    print(model)
#'    print(abs(poles(model$sys)))
#'    print(abs(zeroes(model$sys)))
#' }
test_armamod = function(dim = c(1,1), degrees = c(1,1), b0 = NULL, sigma_L = NULL,
                        digits = NULL, bpoles = NULL, bzeroes = NULL, n.trials = 100) {
  # check input parameter "dim"
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (dim[1] <= 0) || (dim[2] < 0)) {
    stop('argument "dim" must be a vector of integers with length 2, dim[1] > 0 and dim[2] >= 0!')
  }
  # check input parameter "degrees"
  degrees = as.integer(degrees) # note: as.integer converts to vector!
  if ((length(degrees) != 2) || (degrees[1] < 0) || (degrees[2] < (-1))) {
    stop('argument "degrees" must be a vector of integers with length 2, degrees[1] >= 0 and degrees[2] >= -1!')
  }

  m = dim[1]
  n = dim[2]
  p = degrees[1]
  q = degrees[2]

  if (p == 0) bpoles = NULL
  if ( (m != n) || (q <= 0) ) bzeroes = NULL

  # check input parameter "b0"
  if (!is.null(b0)) {
    if ( (!is.numeric(b0)) || (!is.matrix(b0)) || any(dim(b0) != dim) ) {
      stop('"b0" must be a compatible, numeric matrix')
      i = is.na(b0)
      theta = stats::rnorm(sum(i))
      if (!is.null(digits)) theta = round(theta, digits)
      b0[i] = theta
    }
  }
  else {
    b0 = diag(x = 1, nrow = m, ncol = n)
  }

  # check input parameter "sigma_L"
  if (!is.null(sigma_L)) {
    if ( (!is.numeric(sigma_L)) || (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
      stop('"sigma_L" must be a compatible, numeric matrix')
    }
  }
  else {
    sigma_L = matrix(stats::rnorm(n*n), nrow = n, ncol = n)
    if (n >0) {
      sigma_L = t(chol(sigma_L %*% t(sigma_L)))
    }
    if (!is.null(digits)) sigma_L = round(sigma_L, digits)
  }

  i.trial = 0
  err = TRUE
  sd = 1
  while ( (i.trial < n.trials) && (err) ) {
    a = cbind(diag(m), matrix(stats::rnorm(m*m*p, sd = sd), nrow = m, ncol = m*p))
    dim(a) = c(m,m,p+1)
    if (q >= 0) {
      b = cbind(b0, matrix(stats::rnorm(m*n*q, sd = sd), nrow = m, ncol = n*q))
      dim(b) = c(m, n, q+1)
    } else {
      b = array(0, dim = c(m, n, q+1))
    }
    if (!is.null(digits)) {
      a = round(a, digits)
      b = round(b, digits)
    }
    sys = lmfd(a, b)

    err = FALSE
    if ( !is.null(bpoles) ) {
      err = try(min(abs(poles(sys, print_message = FALSE))) <= bpoles, silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    if ( (!err) && (!is.null(bzeroes)) ) {
      err = try((min(abs(zeroes(sys, print_message = FALSE))) <= bzeroes), silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    i.trial = i.trial + 1
    sd = sd/1.1
  }
  if (err) {
    stop('Could not generate a suitable ARMA model with ', n.trials, ' trials!')
  }

  model = armamod(sys = sys, sigma_L = sigma_L)
  return(model)
}
