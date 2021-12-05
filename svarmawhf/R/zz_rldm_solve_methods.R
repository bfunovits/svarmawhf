# Toeplitz Calculations ----

#' Toeplitz Calculations
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' \cr
#' Multiplication of stacked data vector with a block Toeplitz matrix (\code{toepl_fwd()} for MA calculations) or
#' "inversion" of a block Toeplitz matrix in order to perform calculations equivalent to
#' multiplying a given stacked data vector with the inverse of a lower-triangular banded block Toeplitz matrix (\code{toepl_inv()} for AR calculations).
#' Note that matrix polynomials can be mapped one-to-one to banded lower-triangular block Toeplitz matrices.
#'
#' @section MA-type Toeplitz calculations:
#' Given a polynomial matrix of degree \eqn{q} and dimension \eqn{(m \times n)}{(m x n)}, where \eqn{m \geq n}{m \ge n},
#' and given a "wide" input data matrix of dimension \eqn{(n \times nobs)}{(n x nobs)}, where \eqn{nobs} is the number of observations
#' such that each column corresponds to one observation and the number of columns is equal to the number of observations,
#' we calculate a "wide" output data matrix of dimension \eqn{(m \times nobs)}.
#'
#' The function name \code{toepl_fwd} stems from the multiplication of the "stacked" input data vector
#' \eqn{(u_1', \ldots , u_{nobs}')'}{(u[1]', ... , u[nobs]')'} with a banded lower-triangular block Toeplitz matrix \eqn{T}
#' of dimension \eqn{(nobs m \times nobs n)}{(nobs m x nobs n)}
#' whose block elements depend only on the difference between the row- and column-index
#' such that  \eqn{T_{i,j} = d_{i-j}}{T[i,j]=d[i-j]}.
#'
#' @section AR-type Toeplitz calcuations:
#' Given a square polynomial matrix \eqn{c(z)} for which \eqn{c_0}{c[0]} is the identity matrix and given
#' a wide data matrix \code{y = data_in}, obtain the solution \code{u}, a wide data matrix, of the Toeplitz equation
#' \deqn{T (y_1', \ldots , y_{nobs}')' =  (u_1', ... , u_{nobs}')'}{
#'       T (y[1]', \ldots , u[nobs]')' =  (u[1]', ... , u[nobs]')'}
#'
#' Note that the zero-lag coefficient is discarded and the coefficients are in reversed order
#' since this simplifies computations and implementation.
#'
#' @param polm_wide Wide matrix \eqn{( d_0, d_1, \ldots, d_q )}{(d[0], d[1], ... , d[q])}
#'     of dimension \eqn{(m \times (q+1) n)}{(m x (q+1)n)} which represents a
#'     matrix polynomial \eqn{d(z)} of degree \eqn{q}.
#' @param polm_rev Wide matrix \eqn{(c_p, ... , c_1)}{(c[p], ... , c[1])} of dimension \eqn{(n \times (q+1) n)}{(n x p n)}
#'     with coefficients ordered in reverse direction, and
#'     zero-lag coefficient matrix not included.
#'     It represents a square polynomial matrix \eqn{c(z)} with \eqn{c_0}{c[0]} equal to the identity matrix and of degree \eqn{p}.
#' @param data_in Data matrix of dimension \eqn{(dim(inputs) x nobs)}.
#'     Every column corresponds to one observation.
#' @param t0 Integer. Time index where calculations should start.
#'     Default set to 1. For AR calculations, \eqn{degree + 1} would be another smart option.
#'
#' @return data_out Data matrix of dimension \eqn{(dim_out x n_obs)}
#' @export
#' @examples
#' p = test_polm(dim = c(2,2), degree = 2) %>% unclass()
#' dim(p) = c(2,2*3)
#' data = matrix(rnorm(2*100), 2, 100)
#' toepl_fwd(p, data)
toepl_fwd = function(polm_wide, data_in, t0 = 1){
  
  
  # (d_0, d_1, ... , d_q)
  
  # Integer-valued parameters
  dim_out = nrow(polm_wide)
  dim_in = nrow(data_in)
  n_obs = ncol(data_in)
  deg_p = ncol(polm_wide)/dim_in - 1
  
  # Allocate data_out
  data_out = matrix(0, nrow = dim_out, ncol = n_obs)
  
  # Append starting values if t0 is not big enough
  add_zeros = max(deg_p-t0+1, 0)
  if (add_zeros>0){
    data_in = cbind(matrix(0, dim_in, add_zeros),
                    data_in)
    n_obs = n_obs + add_zeros
    t0 = t0 + add_zeros
  }
  
  for (ix in 0:deg_p){
    data_out = data_out + polm_wide[ , (1+ix*dim_in):((ix+1)*dim_in), drop = FALSE] %*% data_in[, (t0-ix):(n_obs-ix), drop = FALSE]
  }
  
  return(data_out)
}


#' @aliases toeplitz_calculations
#' @rdname toepl_fwd
#' @export
#' @examples
#' p = test_polm(dim = c(2,2), degree = 2) %>% unclass()
#' dim(p) = c(2,2*3)
#' data = matrix(rnorm(2*100), 2, 100)
#' toepl_inv(p, data)
toepl_inv = function(polm_rev, data_in, t0 = 1){
  
  # polm = (c_p, ... , c_1)
  
  # Integer-valued parameters
  dim_out = nrow(polm_rev)
  dim_in = nrow(data_in)
  n_obs = ncol(data_in)
  deg_p = ncol(polm_rev)/dim_in
  
  # Fail or return early
  stopifnot(dim_out == dim_in)
  if (deg_p == 0){
    return(data_in)
  }
  
  # Append starting values if t0 is not big enough
  add_zeros = deg_p-t0+1
  if (add_zeros > 0){
    data_in = cbind(matrix(0, dim_in, add_zeros),
                    data_in)
    n_obs = n_obs + add_zeros
    t0 = t0 + add_zeros
  }
  
  if(deg_p>0){
    for (ix_t in t0:n_obs){
      data_in[, ix_t] = data_in[, ix_t, drop = FALSE] - polm_rev %*% matrix(c(data_in[, (ix_t-deg_p):(ix_t-1)]), ncol = 1)
    }
  }
  
  return(data_in[, t0:n_obs, drop = FALSE])
}

# solve_methods.R ----

#' Solve ARMA system
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' \cr
#' Compute the outputs of ARMA(p, q) systems of the form
#' \deqn{y_t = a_1 y_{t-1} + ... + a_p y_{[t-p} + b_0 u_t + \cdots + b_q u_{t-q}}{
#'       y[t] = a[1] y[t-1] + ... + a[p] y[t-p] + b[0] u[t] + ... + b[q] u[t-q]}
#'
#' Values \eqn{y_t}{y[n]}, \eqn{u_t}{u[t]} for \eqn{t\leq 0}{t\le 0} are implicitely set to be zero.
#' However, if we start the iteration with some \eqn{t_0>1}{t0>1} we can enforce non-zero
#' initial values.
#'
#' The routines are used internally and hence do \bold{not} check their arguments.
#' We require \eqn{m > 0}, \eqn{p \geq 0}{p \ge 0}, \eqn{n \geq 0}{n \ge 0}, \eqn{(q+1) \geq 0}{(q+1) \ge 0}
#' and \eqn{1 \leq t_0 \leq N}{1 \le t0 \le N}.
#' Note also that the RcppArmadillo implementation \bold{overwrites} the input argument \code{y}.
#' Use this procedure with care!
#'
#' Note the non standard arguments: The order of the AR coefficients is reversed. The data matrices
#' are organized columnwise (to avoid memory shuffling)!
#'
#' @param a \eqn{(m, mp)} matrix \eqn{(a_p,...,a_1)}{(a[p],...,a[1])}.
#' @param b \eqn{(m, n(q+1))} matrix \eqn{(b_0,...,b_q}{(b[0],...,b[q])}.
#' @param u \eqn{(n, N)} matrix with the inputs \eqn{(u_1,...,u_N}{(u[1],...,u[N])}.
#' @param y \eqn{(m, N)} matrix with the outputs \eqn{(y_1,...,y_N}{(y[1],...,y[N])}.
#' @param t0 integer, start iteration at t=t0.
#'
#' @return The R implementation \code{solve_ARMA_R} returns the matrix \code{y} with the computed outputs.
#'         The RcppArmadillo implementation returns \code{NULL} but \bold{overwrites} the input argument \code{y}!
#'
#' @export
#' @name solve_ARMA
#' @keywords internal
solve_ARMA_R = function(a, b, u, y, t0) {
  m = nrow(y)
  n.obs = ncol(y)
  n = nrow(u)
  p = ncol(a)/m
  q = ncol(b)/n - 1
  if (p > 0) {
    dim(a) = c(m,m,p)
    a = a[,,p:1,drop = FALSE]
  }
  dim(b) = c(m,n,q+1)

  for (t in (t0:n.obs)) {
    y[, t] = matrix(b[,,1], nrow = m, ncol = n) %*% u[,t]
    if (q > 0) {
      for (i in (1:q)) {
        if ((t-i) > 0) y[, t] = y[, t] + matrix(b[ , , i+1], nrow = m, ncol = n) %*% u[, t-i]
      }
    }
    if (p>0) {
      for (i in (1:p)) {
        if ((t-i) > 0) y[, t] = y[, t] + matrix(a[ , , i], nrow = m, ncol = m) %*% y[, t-i]
      }
    }
  }
  return(y)
}


#' Solve (linear) Difference Equations
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' \cr
#' The procedure \code{solve_de} solves the difference equations associated to (V)ARMA models
#' \deqn{a_0 y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p} = b_0 u_t  + b_1 u_{t-1} + ... b_1 u_{t-q}}{
#'       a[0] y[t] + a[1] y[t-1] + ... + a[p] y[t-p] = b[0] u[t]  + b[1] u[t-1] + ... b[q] u[t-q]}
#' or statespace models
#' \deqn{a_{t+1} = A a_t + B u_t \mbox{ and } y_t = C a_t + D u_t.}{
#'       a[t+1] = A a[t] + B u[t] and y[t] = C a[t] + D u[t].}
#'
#' \code{solve_de()} computes the outputs \eqn{y_t}{y[t]} for \eqn{t=1,\ldots,N}{t=1,...,N} for
#' given disturbances \eqn{u_t}{u[t]} \eqn{t=1,\ldots,N}{t=1,...,N}.
#' The starting values  (\eqn{u_t}{u[t]} and \eqn{y_t}{y[t]} for \eqn{t\leq 0}{t \le 0} for VARMA models
#' and \eqn{a_1}{a[1]} for state space models) may be given as optional arguments.
#' The default is to use zero vectors.
#'
#' For the reverse direction, i.e. to reconstruct the disturbances if the outputs are given,
#' the function \code{solve_inverse_de} may be used. In this case the system must be square and
#' the matrix \eqn{D} respectively \eqn{b_0}{b[0]} must be invertible.
#'
#' These functions are mainly intended for internal use and hence only some basic checks
#' on the input parameters are performed.
#'
#' @param sys \code{\link{lmfd}} or \code{\link{stsp}} object
#'            which describes the difference equation.
#' @param u \eqn{(N,n)} matrix with the noise (\eqn{u_t}{u[t]}, \eqn{t=1,...,N}).
#' @param y \eqn{(N,m)} matrix with the outputs (\eqn{y_t}{y[t]}, \eqn{t=1,...,N}).
#' @param u0 \eqn{(h,n)} dimensional matrix with starting values
#'           for the disturbances \eqn{(u_{1-h}, \ldots, u_{-1}, u_0)}{(u[1-h], ..., u[-1], u[0])}.
#'           Note that the last row corresponds to \eqn{u_0}{u[0]}, the last but one row
#'           to \eqn{u_{-1}}{u[-1]} and so on. If \eqn{h>q} then only the last \eqn{q} rows of
#'           \code{u0} are used. In the case \eqn{h<q} the "missing" initial values are set
#'           to zero vectors. \cr
#'           The default value \code{u0=NULL} sets all initial values \eqn{u_t}{u[t]}, \eqn{t \leq 0}{t \le 0}
#'           equal to zero vectors.
#' @param y0 \eqn{(h,m)} dimensional matrix with starting values
#'           for the outputs \eqn{(y_{1-h}, \ldots, y_{-1}, y_0)}{(y[1-h], ..., y[-1], y[0])}.
#'           This (optional) parameter is interpreted analogously to \code{u0}.
#' @param ... not used.
#'
#'
#' @return List with slots
#' \item{y}{\eqn{(N,n)} matrix with the (computed) outputs.}
#' \item{u}{\eqn{(N,n)} matrix with the (computed) noise.}
#' \item{a}{\eqn{(N+1,n)} matrix with the (computed) states (\eqn{a_t}{a[t]}, \eqn{t=1,...,N+1}).
#'          Note that this matrix has (\eqn{N+1}) rows! This slot is only present for state space models.}
#' @export
#'
#' @examples
#' ### generate a random ARMA(2,1) model (with two outputs) #########
#' model = test_armamod(dim = c(2,2), degrees = c(2,1),
#'                      digits = 2, bpoles = 1, bzeroes = 1)
#'
#' # generate random noise sequence (sample size N = 100)
#' u = matrix(rnorm(100*2), nrow = 100, ncol = 2)
#'
#' # generate random initial values
#' u0 = matrix(rnorm(2), nrow = 1, ncol = 2) # u[0]
#' y0 = matrix(rnorm(2), nrow = 1, ncol = 2) # y[0]
#'
#' # compute outputs "y[t]"
#' # note that y0 has only one row, thus y[-1] is set to zero!
#' data = solve_de(model$sys, u = u, y0 = y0, u0 = u0)
#'
#' # we can reconstruct the noise "u" from given outputs "y"
#' data1 = solve_inverse_de(model$sys, y = data$y, u0 = u0, y0 = y0)
#' all.equal(data$u, data1$u)
#'
solve_de = function(sys, u, ...)  {
  UseMethod("solve_de", sys)
}


#' @rdname solve_de
#' @export
solve_de.lmfd = function(sys, u, u0 = NULL, y0 = NULL, ...) {

  a = unclass(sys$a)
  b = unclass(sys$b)
  m = dim(a)[1]
  p = dim(a)[3] - 1
  n = dim(b)[2]
  q = dim(b)[3] - 1

  if ( (m*(p+1)) == 0 ) stop('illegal ARMA system (m=0 or p<0)')

  # check 'u'
  if ( (!is.matrix(u)) || (ncol(u) != n) ) stop('disturbances "u" are not compatible to the model!')
  n.obs = nrow(u)

  # check starting values u0
  if (is.null(u0)) {
    u0 = matrix(0, nrow = max(p,q), ncol = n)
  }
  if ( (!is.matrix(u0)) || (ncol(u0) != n) ) stop('initial disturbances "u0" are not compatible to the model!')
  # just keep the last max(p,q) rows / prepend zero rows if necessary
  qq = min(nrow(u0), max(p,q))
  u0 = rbind(matrix(0, nrow = max(p,q) - qq, ncol = n), u0[iseq(nrow(u0)-qq+1,nrow(u0)),,drop = FALSE])

  # check starting values y0
  if (is.null(y0)) {
    y0 = matrix(0, nrow = max(p,q), ncol = m)
  }
  if ( (!is.matrix(y0)) || (ncol(y0) != m) ) stop('initial outputs "y0" are not compatible to the model!')
  # just keep the last max(p,q) rows / prepend zero rows if necessary
  pp = min(nrow(y0), max(p,q))
  y0 = rbind(matrix(0, nrow = max(p,q) - pp, ncol = m), y0[iseq(nrow(y0)-pp+1, nrow(y0)),,drop = FALSE])

  if (n.obs == 0) {
    return(list(y = matrix(0, nrow = n.obs, ncol = m), u = u))
  }

  # convert ARMA system
  #    a[0] y[t] + a[1] y[t-1] + ... = b[0] u[t] + b[1] u[t-1] + ...
  # to system of the form
  #    y[t] = a[1] y[t-1] + ... + b[0] u[t] ...
  # and create parameter matrices a = (a[p], ..., a[1]) and b = (b[0], ..., b[q])
  a0 = matrix(a[,,1], nrow = m, ncol = m)

  dim(b) = c(m,n*(q+1))
  if ((n*(q+1)) > 0) {
    b = solve(a0, b)
  }

  # note for solve_ARMA_cpp we have to reshuffle the AR parameter as:  a = (a[p],...,a[1])
  if (p>0) {
    a = a[,,(p+1):2, drop = FALSE]
    dim(a) = c(m, m*p)
    a = -solve(a0, a)
  } else {
    a = matrix(0, nrow = m, ncol = 0)
  }


  # transpose u matrix and prepend initial values
  u = t(rbind(u0,u))     # (n, max(p,q)+n.obs)

  # create y matrix (including starting values)
  y = t(rbind(y0, matrix(0, nrow = n.obs, ncol = m)))  # (m, max(p,q)+n.obs)

  # solve ARMA system
  .Call(`_svarmawhf_solve_ARMA_cpp`, a, b, u, y, max(p,q)+1)

  y = t(y[, (max(p,q)+1):(n.obs+max(p,q)), drop = FALSE])
  u = t(u[, (max(p,q)+1):(n.obs+max(p,q)), drop = FALSE])

  return(list(y = y, u = u))
}


#' @rdname solve_de
#' @export
solve_inverse_de = function(sys, y, ...)  {
  UseMethod("solve_inverse_de", sys)
}

#' @rdname solve_de
#' @export
solve_inverse_de.lmfd = function(sys, y, u0 = NULL, y0 = NULL, ...) {

  a = unclass(sys$a)
  b = unclass(sys$b)
  m = dim(a)[1]
  p = dim(a)[3] - 1
  n = dim(b)[2]
  q = dim(b)[3] - 1

  if ( (m != n) || (min(m, p+1, q+1) <= 0) ) {
    stop('"solve_inverse_de" is only implemented for square systems (m = n > 0 and p,q >=0)!')
  }

  # check 'y'
  if ( (!is.matrix(y)) || (ncol(y) != m) ) stop('outputs "y" are not compatible to the model!')
  n.obs = nrow(y)

  # check starting values u0
  if (is.null(u0)) {
    u0 = matrix(0, nrow = max(p,q), ncol = n)
  }
  if ( (!is.matrix(u0)) || (ncol(u0) != n) ) stop('initial disturbances "u0" are not compatible to the model!')
  # just keep the last max(p,q) rows / prepend zero rows if necessary
  qq = min(nrow(u0), max(p,q))
  u0 = rbind(matrix(0, nrow = max(p,q) - qq, ncol = n), u0[iseq(nrow(u0)-qq+1,nrow(u0)),,drop = FALSE])

  # check starting values y0
  if (is.null(y0)) {
    y0 = matrix(0, nrow = max(p,q), ncol = m)
  }
  if ( (!is.matrix(y0)) || (ncol(y0) != m) ) stop('initial outputs "y0" are not compatible to the model!')
  # just keep the last max(p,q) rows / prepend zero rows if necessary
  pp = min(nrow(y0), max(p,q))
  y0 = rbind(matrix(0, nrow = max(p,q) - pp, ncol = m), y0[iseq(nrow(y0)-pp+1, nrow(y0)),,drop = FALSE])

  if (n.obs == 0) {
    return(list(y = y, u = matrix(0, nrow = n.obs, ncol = n)))
  }

  # convert ARMA system
  #    a[0] y[t] + a[1] y[t-1] + ... = b[0] u[t] + b[1] u[t-1] + ...
  # to system of the form
  #    u[t] = b[1] u[t-1] + ... + a[0] y[t] ...
  # and create parameter matrices b = (b[q], ..., b[1]) and a = (a[0], ..., a[p])
  b0 = matrix(b[,,1], nrow = m, ncol = m)

  dim(a) = c(m,m*(p+1))
  if ((m*(p+1)) > 0) {
    a = solve(b0, a)
  }

  # note for solve_ARMA_cpp we have to reshuffle the MA parameter as:  b = (b[q],...,q[1])
  if (q>0) {
    b = b[,,(q+1):2, drop = FALSE]
    dim(b) = c(m, m*q)
    b = -solve(b0, b)
  } else {
    b = matrix(0, nrow = m, ncol = 0)
  }

  # transpose y matrix and prepend initial values
  y = t(rbind(y0, y))     # (m, max(p,q)+n.obs)

  # create u matrix (including starting values)
  u = t(rbind(u0, matrix(0, nrow = n.obs, ncol = n)))  # (n, max(p,q)+n.obs)

  # solve ARMA system
  .Call(`_svarmawhf_solve_ARMA_cpp`, b, a, y, u, max(p,q)+1)

  y = t(y[, (max(p,q)+1):(n.obs+max(p,q)), drop = FALSE])
  u = t(u[, (max(p,q)+1):(n.obs+max(p,q)), drop = FALSE])

  return(list(y = y, u = u))
}
