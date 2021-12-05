#' Autocovariance, Autocorelation and Partial Autocorrelation Function
#'
#' This function was originally part of R package \strong{RLDM}.
#' \cr
#' Compute respectively estimate the autocovariance, autocorrelation or partial autocorrelation function of a stationary process.
#'
#' The class of the input parameter "\code{obj}" determines the S3 method called and hence what is actually computed.
#'
#' \bold{Population ACF:}
#'
#' If "\code{obj}" is an \code{\link{armamod}} or \code{\link{stspmod}} object then \code{autocov(obj, ...)} computes the ACF of the corresponding stationary process.
#'
#' Note however, that the function returns nonsense, if the model does not satisfy the \emph{stability} condition.
#'
#' \bold{Change the type of an ACF:}
#'
#' Calling \code{autocov(obj, type)}, where "\code{obj}" is an \code{autocov} object returns an ACF of the desired type.
#' E.g. if "\code{obj}" holds a partial autocorrelation function then \code{autocov(obj, type = 'covariance')} may be used to retrieve the corresponding autocovariance function.
#'
#' This is possible since the \code{autocov} object stores the "original" autocovariances in a slot named \code{gamma}.
#'
#' \bold{Sample ACF:}
#'
#' The default S3 method estimates the ACF from given data.
#' It assumes that "\code{obj}" is a univariate or multivariate numeric time series object,
#' a (numeric) data frame or a (numeric) vector, respectively matrix
#' and then simply calls the function \code{\link[stats]{acf}} in the \pkg{stats} package to compute the sample autocovariance function.
#' If needed, then the corresponding sample autocorrelation, respectively sample partial autocorrelation function is computed (and returned).
#'
#' The syntax is quite analogous to \code{\link[stats]{acf}}, so please consider the documentation of \code{\link[stats]{acf}} for more details.
#'
#' Note that \pkg{stats} stores autocovariance/autocorrelation functions as \code{(lag.max+1,m,m)} dimensional arrays, whereas \pkg{RLDM} uses \code{(m,m,lag.max+1)} dimensional arrays.
#'
#' The definition of partial autocorrelations used by \code{\link[stats]{acf}} differs from the definition used here.
#' Furthermore \code{\link[stats]{acf}} skips the lag zero partial autocorrelation coefficient
#' and thus the pacf computed by  \code{\link[stats]{acf}} is (lag.max,n,n) dimensional.
#'
#' The default choice for the number of lags is \eqn{10*log10(N/m)} where \eqn{N}
#' is the number of observations and \eqn{m} the number of series.
#' This number will be automatically limited to one less than the number
#' of observations in the series.
#'
#' @param obj either a \code{\link{armamod}}, \code{\link{stspmod}}, \code{\link{autocov}} object
#'            or a "data" object.
#' @param type   character string giving the type of acf to be computed. Allowed values are
#'               "covariance" (the default), "correlation", or "partial". Will be partially matched.
#'               Note that the default value here is "covariance" whereas \code{\link[stats]{acf}}
#'               uses "correlation" as default.
#' @param ...    not used.
#' @param lag.max (integer) maximum lag.
#' @param na.action function to be called to handle missing values. \code{\link[stats]{na.pass}} can be used.
#' @param demean logical. Should the covariances be about the sample means?
#'
#' @return \code{autocov} object, i.e. a list with slots
#' \item{acf}{\code{\link[rationalmatrices]{pseries}} object, which stores the covariances (correlations).}
#' \item{type}{character string which indicates the type of the ACF.}
#' \item{gamma}{(m,m,lag.max+1) dimensional array which stores the autocovariance function.}
#' \item{names}{(m)-dimensional character vector or NULL. This optional slot stores the names
#'              for the components of the time series/process.}
#' \item{label}{character string or NULL.}
#' \item{n.obs}{integer or NULL. This slot stores the sample size.}
#'
#'
#' @seealso The autocovariance function of (V)ARMA processes may also be
#' computed by \code{\link[stats]{ARMAacf}} in the scalar case and by
#' \code{\link[MTS]{VARMAcov}} in the multivariate case (m > 1).
#'
#' As noted above the sample ACF is computed via the \code{\link[stats]{acf}} routine in the \pkg{stats} package.
#'
#' @export
autocov = function(obj, type, ...) {
  UseMethod("autocov", obj)
}

#' @rdname autocov
#' @export
autocov.default = function(obj, type=c('covariance','correlation','partial'), lag.max = NULL,
                           na.action = stats::na.fail, demean = TRUE, ...) {
  
  if (!is.null(lag.max)) {
    lag.max = as.integer(lag.max)[1]
    if (lag.max < 0) stop('negative lag.max')
  }
  type = match.arg(type)
  
  out_acf = try(stats::acf(obj, lag.max = lag.max, type = 'covariance', plot = FALSE,
                           na.action = na.action, demean = demean))
  
  if (inherits(out_acf, 'try-error')) stop('stats:acf failed, the input "obj" may not be supported.')
  
  gamma = aperm(out_acf$acf, c(2,3,1))
  n.obs = out_acf$n.used
  
  if (type == 'covariance') {
    acf = gamma
  }
  if (type == 'correlation') {
    # compute correlations
    m = dim(gamma)[1]
    # standard deviations
    s = sqrt(diag(matrix(gamma[,,1], nrow = m, ncol = m)))
    s = matrix(s, nrow = m, ncol = m)
    s2 = s * t(s)
    acf = gamma
    for (i in (1: dim(gamma)[3])) acf[,,i] = acf[,,i] / s2
  }
  if (type == 'partial') {
    acf = est_ar_dlw(gamma, penalty = -1)$partial
  }
  acf = structure(acf, class = c('pseries', 'ratm') )
  
  out = structure( list(acf = acf, type = type, gamma = gamma, names = NULL, label = NULL,  n.obs = n.obs),
                   class = c('autocov', 'rldm') )
  
  return(out)
}

