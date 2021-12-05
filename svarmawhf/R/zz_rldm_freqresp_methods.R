# freqresp_methods.R ######################################
# defines the main methods for the 'freqresp' class

#' Discrete Time Fourier Transform
#'
#' Compute the Discrete Time Fourier Transform for data
#' stored in a 3-dimensional array.
#'
#' @param a \eqn{(m,n,k)} dimensional (numeric) array
#' @param n.f (integer) number of frequencies
#'
#' @return A \code{\link[rationalmatrices]{zvalues}} object.
#' @export
#'
#' @keywords internal
dft_3D = function(a, n.f = dim(a)[3]) {
  n.f = as.integer(n.f)[1]
  if (n.f > 0) {
    z = exp((0:(n.f-1)) * (2*pi*complex(imaginary = -1)/n.f))
  } else {
    z = complex(0)
  }

  d = dim(a)
  if (min(c(d,n.f)) == 0) {
    a = array(complex(real = 0), dim = c(d[1], d[2], n.f))
    attr(a,'z') = z
    class(a) = c('zvalues','ratm')
    return(a)
  }

  dim(a) = c(d[1]*d[2], d[3])
  a = t(a)
  if (d[3] > n.f) {
    h = ceiling(d[3]/n.f)
    a = rbind(a, matrix(0, nrow = h*n.f - d[3], ncol = d[1]*d[2]))
    dim(a) = c(n.f, h, d[1]*d[2])
    a = apply(a, MARGIN = c(1,3), FUN = sum)
  }
  if (d[3] < n.f) {
    a = rbind(a, matrix(0, nrow = n.f - d[3], ncol = d[1]*d[2]))
  }
  # print(dim(a))
  a = stats::mvfft(a)
  # print(dim(a))
  a = t(a)
  dim(a) = c(d[1],d[2],n.f)
  attr(a,'z') = exp((0:(n.f-1)) * (2*pi*complex(imaginary = -1)/n.f))
  class(a) = c('zvalues','ratm')
  return(a)
}


#' Frequency Response Function
#'
#' Compute the \emph{frequency response function} (also called \emph{transfer function}) associated
#' to a VARMA or statespace model.
#'
#' The frequency response function (or transfer function) associated
#' to an ARMA or statespace model is
#' \deqn{
#' K(\lambda) = \sum_{j=0}^{\infty} k_j e^{-i\lambda j}
#' }{
#' K(\lambda) = sum_{j=0}^{\infty} k[j] exp(-i\lambda j)
#' }
#' where \eqn{(k_j \,|\, j\geq 0)}{(k[j], j \ge 0)} is the \emph{impulse response} of the model.
#' See also \code{\link{impresp}}.
#'
#' For an ARMA model the frequency response is equal to
#' \deqn{
#' K(\lambda) =
#' (a_0 + a_1 e^{-i\lambda} + \cdots + a_p e^{-i\lambda p})^{-1}
#' (b_0 + b_1 e^{-i\lambda} + \cdots + b_q e^{-i\lambda q})
#' }{
#' K(\lambda) =
#' (a[0] + a[1] exp(-i\lambda) + \dots + a[p] exp(-i\lambda p))^{-1}
#' (b[0] + b[1] exp(-i\lambda) + \dots + b[q] exp(-i\lambda q))
#' }
#' and for a statespace model we have
#' \deqn{
#' K(\lambda) = C(e^{i\lambda}I_s - A)^{-1}B+D
#' }{
#' K(\lambda) = C(exp(i\lambda) I - A)^{-1} B + D
#' }
#' Note that \eqn{K()} is the discrete-time Fourier transform (DTFT) of the impulse response.
#' If the impulse response is absolutely summable then the coefficents \eqn{k_j}{k[j]}
#' may be reconstructed from the frequency response via the inverse DTFT
#' \deqn{
#' k_j = \frac{1}{2\pi} \int_{-\pi}^{\pi} K(\lambda) e^{i\lambda j} d\lambda
#' }{
#' k[j] = (1/2\pi) int_{-\pi}^{\pi} K(\lambda) e^{i\lambda j} d\lambda
#' }
#'
#' The S3 methods \code{freqresp.*} evaluate the function on a grid of angular
#' frequencies \eqn{\lambda_j = 2\pi j/N}{\lambda[j] = 2\pi j/N}, \eqn{j=0,\ldots,N-1}{j=0,\dots,N-1}
#' and store the result (together with \code{sigma_L}) in a \strong{freqresp} object.
#'
#'
#' @param obj \code{\link{armamod}}, \code{\link{stspmod}} or \code{\link{impresp}} object. Note that
#'            for an impulse response object the result is only an approximation
#'            of the "true" frequency response due to the finite number of coefficients.
#' @param n.f number of frequencies.
#' @param ... not used.
#'
#' @return \code{freqresp} object, i.e. a list with slots
#' \item{frr}{\code{\link[rationalmatrices]{zvalues}} object.}
#' \item{sigma_L}{(n,n)-dimensional matrix which contains a left square root
#'                of the noise covariance matrix \eqn{\Sigma}.}
#' \item{names}{(m)-dimensional character vector or NULL. This optional slot stores the names
#'              for the components of the time series/process.}
#' \item{label}{character string or NULL.}
#'
#' @export
freqresp = function(obj, n.f, ...) {
  UseMethod("freqresp", obj)
}

#' @rdname freqresp
#' @export
freqresp.armamod = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  d = dim(obj$sys)
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  if (n.f > 0) {
    z = exp((0:(n.f-1)) * (2*pi*complex(imaginary = -1)/n.f))
  } else {
    z = complex(0)
  }


  if ( min(c(m,n,q+1,n.f)) == 0) {
    frr = array(complex(real = 0), dim = c(m, n, n.f))
    attr(frr, 'z') = complex(0)
    class(frr) = c('zvalues', 'ratm')

    out = structure(list(frr = frr, sigma_L = obj$sigma_L, names = obj$names, label = obj$label),
                    class = c('freqresp', 'rldm'))
  }
  a = unclass(dft_3D(unclass(obj$sys$a), n.f))
  frr = unclass(dft_3D(unclass(obj$sys$b), n.f))
  if (m == 1) {
    a = aperm(array(a, dim = c(m, n.f, n)), c(1, 3, 2))
    frr = frr / a
  } else {
    for (i in (1:n.f)) {
      b = try(solve(a[,,i], frr[,,i]), silent = TRUE)
      if (!inherits(b, 'try-error')) {
        frr[,,i] = b
      } else {
        frr[,,i] = NA_complex_
      }
    }
  }
  attr(frr, 'z') = z
  class(frr) = c('zvalues', 'ratm')

  out = structure(list(frr = frr, sigma_L = obj$sigma_L, names = obj$names, label = obj$label),
                  class = c('freqresp', 'rldm'))
  return(out)
}


#' @rdname freqresp
#' @export
freqresp.stspmod = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  frr = zvalues(obj$sys, n.f = n.f)
  # make sure that the frr object is of type 'complex'
  if (!is.complex(frr)) {
    frr = unclass(frr)
    z = attr(frr, 'z')
    d = dim(frr)
    frr = array(as.complex(frr), dim = d)
    frr = structure(frr, z = z, class = c('zvalues', 'ratm'))
  }

  out = structure(list(frr = frr, sigma_L = obj$sigma_L, names = obj$names, label = obj$label),
                  class = c('freqresp', 'rldm'))
  return(out)
}




#' @rdname freqresp
#' @export
freqresp.impresp = function(obj, n.f = 128, ...) {
  n.f = as.integer(n.f)[1]
  if (n.f < 0) stop('the number of frequencies "n.f" must be a non negative integer')

  irf = unclass(obj$irf)
  frr = dft_3D(irf, n.f)

  out = structure(list(frr = frr, sigma_L = obj$sigma_L, names = obj$names, label = obj$label),
                  class = c('freqresp', 'rldm'))
  return(out)
}


