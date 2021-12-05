#' Simulate from a State Space or VARMA Model
#'
#' Generate a time series from a given process model.
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' \cr
#' In order to generate a "stationary" trajectory (of a stable model) one
#' has to chose suitable initial starting values.
#' This is not quite easy, in particular for VARMA models.
#' As a simple remedy, the procedure offers the option of a "burn-in" phase.
#' The length of this phase has to be chosen by the user.
#'
#' If the user would like to have more control on the disturbances
#' and the initial values, then \code{\link{solve_de}} may be used.
#'
#' @param model either a \code{\link{armamod}} or \code{\link{stspmod}} object.
#' @param n.obs sample size (\eqn{N}).
#' @param rand.gen an (optional) function to generate the disturbances \eqn{u_t}{u[t]}. Note that
#'              \code{rand.gen()} should generate an iid sample of a random variable with mean zero
#'              and variance one.
#' @param n.burnin length of an initial "burn-in" phase (denoted with \eqn{N_0}{N0}).
#' @param ...   not used.
#'
#' @return List with slots
#' \item{y}{\eqn{(N,m)} matrix with the generated outputs.}
#' \item{u}{\eqn{(N,n)} matrix with the noise.}
#' \item{a}{\eqn{(N+1,s)} matrix with the generated states (\eqn{a_t}{a[t]}, \eqn{t=1,...,N+1}).
#'          Note that this matrix has (\eqn{N+1}) rows! This slot is only present for state space models.}
#'
#' @export
#'
#' @examples
#' # Random Walk ############################################################################
#' model = armamod(lmfd(a = c(1,-1), b = 1))
#' # generate outputs "y"
#' n.obs = 100
#' data = sim(model, n.obs = n.obs, y0 = 1)
#' plot(data$y, type = 'l')
#'
#' # bivariate ARMA(2,1) model ##############################################################
#' model = test_armamod(dim = c(2,2), degrees = c(2,1), bpoles = 1, bzeroes = 1)
#' # generate outputs "y" with zero initial conditions
#' n.obs = 50
#' data =  sim(model, n.obs = n.obs)
#' # reconstruct noise "u" from given outputs "y"
#' data1 = solve_inverse_de(model$sys, y = data$y)
#' all.equal(data$u, data1$u)
#'
sim = function(model, n.obs, rand.gen, n.burnin, ...) {
  UseMethod("sim", model)
}

#' @rdname sim
#' @export
sim.armamod = function(model, n.obs, rand.gen = stats::rnorm, n.burnin = 0, ...) {

  d = dim(model$sys)
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  if ( (m*(p+1)) == 0 ) stop('illegal ARMA system (m=0 or p<0)')

  # check 'n.obs'
  n.obs = as.integer(n.obs[1])
  if (n.obs <= 0) stop('"n.obs" must be a positive integer!')

  # check 'n.burnin'
  n.burnin = as.integer(n.burnin[1])
  if (n.burnin < 0) stop('"n.burnin" must be a non negative integer!')

  # get ARMA parameters
  a = unclass(model$sys$a)
  b = unclass(model$sys$b)

  # convert ARMA system
  #    a[0] y[t] + a[1] y[t-1] + ... = b[0] u[t] + b[1] u[t-1] + ...
  # to system of the form
  #    y[t] = a[1] y[t-1] + ... + b[0] u[t] ...
  # and create parameter matrices a = (a[p], ..., a[1]) and b = (b[0], ..., b[q])
  a0 = matrix(a[,,1], nrow = m, ncol = m)

  dim(b) = c(m, n*(q+1))
  if ((n*(q+1)) > 0) {
    b = solve(a0, b)
  }

  # note for solve_ARMA_cpp we have to reshuffle the AR parameter as:  a = (a[p],...,a[1])
  if (p > 0) {
    a = a[,,(p+1):2, drop = FALSE]
    dim(a) = c(m, m*p)
    a = -solve(a0, a)
  } else {
    a = matrix(0, nrow = m, ncol = 0)
  }

  # generate disturbances
  u = model$sigma_L %*% matrix(rand.gen((n.obs+n.burnin)*n), nrow = n, ncol = n.obs+n.burnin)
  # outputs
  y = matrix(0, nrow = m, ncol = n.obs+n.burnin)

  # solve ARMA system
  .Call(`_svarmawhf_solve_ARMA_cpp`, a, b, u, y, 1)

  if (n.burnin > 0) {
    y = y[, (n.burnin+1):(n.burnin+n.obs), drop = FALSE]
    u = u[, (n.burnin+1):(n.burnin+n.obs), drop = FALSE]
  }

  return(list(y = t(y), u = t(u)))
}

