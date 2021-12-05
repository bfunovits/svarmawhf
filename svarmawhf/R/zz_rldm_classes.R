# classes.R ######################################
# defines the main classes


#' Creator for MFD classes
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' \cr
#' Either a left-matrix fraction description (MFD) or a right-MFD plus parameterisation of noise covariance.
#'
#' In Hannan, Deistler (2012, page 7), RMFDs are also called dynamic adjustment forms.
#' Internally, MFDs are lists with slots \code{sys}, \code{sigma_L}, \code{names}, \code{label}.
#'
#' @param sys \code{\link{lmfd}} or \code{\link{rmfd}} object
#' @param sigma_L Left-factor of noise covariance,
#'   i.e. the covariance \eqn{\sigma} is obtained as \code{sigma_L * t(sigma_L)}.
#'   If \code{sigma_L} is a vector of dimension \eqn{n}, where \eqn{n} is the input dimension, only the diagonal elements are parametrized.
#'   If it is a vector of dimension \eqn{n^2}, then the elements of \code{sigma_L} are filled column by column.
#' @param names optional vector of character strings
#' @param label optional character string
#'
#' @return Object of class \code{armamod}.
#'
#' @export
#'
#' @examples
#' x = armamod(sys = lmfd(c(1, 0.5), 1), sigma_L = diag(1))
#' x
#'
#' @name mfds
armamod = function(sys, sigma_L = NULL, names = NULL, label = NULL) {
  if (!inherits(sys, 'lmfd')) stop('"sys" must be an lmfd object')

  d = dim(sys)
  m = d[1]
  n = d[2]

  if (is.null(sigma_L)) sigma_L = diag(n)
  if (!is.numeric(sigma_L)) stop('parameter sigma_L is not numeric')
  if ( is.vector(sigma_L) ) {
    if (length(sigma_L) == n) sigma_L = diag(sigma_L, nrow = n, ncol = n)
    if (length(sigma_L) == (n^2)) sigma_L = matrix(sigma_L, nrow = n, ncol = n)
  }
  if ( (!is.matrix(sigma_L)) || any(dim(sigma_L) != n) ) {
    stop('"sigma_L" is not compatible')
  }

  x = list(sys = sys, sigma_L = sigma_L, names = names, label = label)
  x = structure(x, class = c('armamod', 'rldm'))
  return(x)
}
