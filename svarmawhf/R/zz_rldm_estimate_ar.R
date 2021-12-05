#' Estimate Autoregressive Models
#'
#' This function was originally part of R package \strong{RLDM}.
#' \cr
#' The function \code{est_ar} estimates (V)AR models
#' \deqn{(y_t - \mu) = a_1 (y_{t-1} - \mu) + \cdots + a_p (y_{t-p} - \mu) + u_t}{
#'        (y[t] - \mu) = a[1] (y[t-1] - \mu) + ... + a[p] (y[t-p] - \mu) + u[t]}
#' from a given sample or a given (sample) autocovariance function.
#' The model order \eqn{p} is chosen by an information criterion, like \emph{AIC} or \emph{BIC}.
#' The "helper" functions \code{est_ar_ols}, \code{est_ar_yw} and \code{est_ar_dlw} implement
#' the three available estimation methods: estimation by ordinary least squares, the Yule-Walker
#' estimates and the Durbin-Levinson-Whittle method.
#'
#' @section OLS method:
#'
#' The helper function \code{est_ar_ols} implements three schemes to estimate the mean \eqn{\mu}
#' and the AR parameters. The choice \code{mean_estimate = "zero"} assumes \eqn{\mu=0} and thus
#' the AR parameters are determined from the regression:
#' \deqn{y_t = a_1 y_{t-1} + \cdots + a_p y_{t-p} + u_t \mbox{ for } t=p+1,\ldots,N}{
#'       y[t] = a[1] y[t-1] + ... + a[p] y[t-p] + u[t] for t=p+1,...,N}
#' In the case \code{mean_estimate = "sample.mean"} the mean \eqn{\mu} is estimated by the sample
#' mean and the AR parameters are determined by the LS estimate of the regression
#' \deqn{(y_t - \mu) = a_1 (y_{t-1} - \mu) + \cdots + a_p (y_{t-p} - \mu) + u_t \mbox{ for } t=p+1,\ldots,N}{
#'       (y[t] - \mu) = a[1] (y[t-1] - \mu) + ... + a[p] (y[t-p] - \mu) + u[t] for t=p+1,...,N}
#' In the last case \code{mean_estimate = "intercept"}, a regression with intercept
#' \deqn{y_t = d + a_1 y_{t-1} + \cdots + a_p y_{t-p} + u_t \mbox{ for } t=p+1,\ldots,N}{
#'       y[t] = d + a[1] y[t-1] + ... + a[p] y[t-p] + u[t] for t=p+1,...,N}
#' is considered. The estimate for \eqn{\mu} then is obtained as
#' \deqn{\mu = (I_m - a_1 - \cdots - a_p)^{-1} d}{
#'       \mu = (I - a[1] - ... - a[p])^{-1} d}
#' This estimate of the mean \eqn{\mu} fails if the estimated AR model has a \emph{unit root}, i.e. if
#' \eqn{(I_m - a_1 - \cdots - a_p)}{(I - a[1] - ... - a[p])} is singular.
#'
#' The sample covariance of the corresponding residuals (scaled with \eqn{1/(N-p)}) serves as
#' an estimate for the noise covariance \eqn{\Sigma}.
#'
#' For the actual computations the routine \code{\link[stats]{lsfit}} in the \pkg{stats} package is used.
#'
#' @section Yule-Walker estimates:
#'
#' Both  \code{est_ar_yw} and \code{est_ar_dlw} use the \emph{Yule-Walker} equations to estimate
#' the AR coefficients \eqn{(a_i)}{(a[i])} and the noise covariance matrix \eqn{\Sigma}.
#' However, they use a different numerical scheme to solve these equations.
#' The function \code{est_ar_dlw} uses the \emph{Durbin-Levinson-Whittle} recursions and
#' in addition returns (estimates of) the partial autocorrelation coefficients.
#'
#' If \code{obj} is a "time series" object, then first the ACF is estimated with a call to
#' \code{\link{autocov}}. The option \code{mean_estimate = "zero"} implies that the mean is assumed to
#' be zero (\eqn{\mu = 0}) and therefore \code{autocov} is called with the option \code{demean = FALSE}.
#'
#' For \code{mean_estimate = "sample.mean"} or \code{mean_estimate = "intercept"} the mean \eqn{\mu}
#' is estimated by the sample mean and the ACF is computed with  \code{demean = TRUE}.
#'
#' @section Estimation of the AR order:
#'
#' The order \eqn{p} of the AR model is chosen by minimizing an information criterion of the form
#' \deqn{IC(p) = \ln\det\Sigma_p + c(p)r(N) \mbox{ for } p = 0,\ldots,p_{\max}}{
#'       IC(p) = ln(det(\Sigma[p])) + c(p)r(N) for p = 0,...,p.max}
#' where \eqn{\Sigma_p}{\Sigma[p]} is the estimate of the noise (innovation) covariance,
#' \eqn{c(p)} counts the number of parameters of the model,
#' and \eqn{r(N)} is the "penalty" per parameter of the model.
#' Note that \eqn{\log\det\Sigma}{log(det(\Sigma))} is up to a constant
#' and a scaling factor \eqn{-(N-p)/2}
#' equal to the (scaled, approximate) Gaussian log likelihood of the model
#' \deqn{ll = -(1/2)(m \ln(2\pi) + m + \ln\det \Sigma_p)}{
#'       ll = -(1/2)(m ln(2\pi) + m + ln(det(\Sigma[p])))}
#' See also \code{\link{ll}}. Note that the value \eqn{ll}, which is returned by
#' this routine, is the (approximate) log Likelihood **scaled** by a factor \eqn{1/(N-p)}.
#'
#' For an AR(p) model with intercept, the number of parameters is
#' \eqn{c(p) = p m^2 + m}{c(p) = p*m^2 + m} and for the AR model without intercept
#' \eqn{c(p) = p m^2}{c(p) = p*m^2}.
#'
#' The Akaike information criterion (AIC) corresponds to \eqn{r(N)=2/N} and the
#' Bayes information criterion (BIC) uses the penatlty \eqn{r(N)=\log(N)/N}{r(N)=log(N)/N}.
#'
#' For the helper routines, the user has to set the penalty term \eqn{r(N)} explicitly via the
#' input parameter "\code{penalty}". The default choice \code{penalty = -1} means that the maximum
#' possible order \code{p=p.max} is chosen.
#'
#' The function \code{est_ar} offers the parameter "\code{ic}" which tells the routine to set the
#' penalty accordingly. Note that the choice \code{ic="max"} sets \eqn{r(N) = -1} and thus again
#' the model with maximum possible order is fitted.
#'
#' The default maximum order \code{p.max} is chosen as follows.
#' The helper functions \code{est_ar_yw} and \code{est_ar_dlw} simply chose the maximum
#' accordng to maximum lag of the given autocovariances, \code{p.max = dim(gamma)[3] - 1}.
#' The routine \code{est_ar_ols} uses the minimum of \eqn{12}, \eqn{(N-1)/(m+1)} and \eqn{10*log10(N)}
#' as default. The function \code{est_ar} uses the same value.
#' However, if "\code{obj}" is an \code{autocov} object then \code{p.max} is in addition bounded
#' by the number of lags contained in this object.
#'
#' @section Notes:
#'
#' The Yule-Walker estimates offer an easy way to reconstruct the "true" model if the population
#' autocovariance function is given. The noise covariance (and thus the likelihood values) should
#' not improve when a model with an order larger than the true model order is "estimated".
#' However due to numerical errors this may not be true. As a simple trick one may call \code{est_ar}
#' (\code{est_ar_yw} or \code{est_ar_dlw}) with a very small positive \code{penalty}. See the example below.
#'
#' The functions are essentially equivalent to the \pkg{stats} routines.
#' They are (re) implemented for convenience, such that the input and output parameters (models) fit
#' to the \pkg{RLDM} conventions.
#'
#' The AIC values of \pkg{RLDM} routines are equivalent to the AIC values computed
#' by the \pkg{stats} routines up to a constant and up to scaling by \eqn{N}.
#'
#' It seems that the Yule-Walker estimate \code{stats::\link[stats]{ar.yw}} uses a scaling
#' factor \eqn{(N - m(p+1))/N} for the noise covariance \eqn{\Sigma}.
#'
#' Finally note that \code{est_ar_ols}, \code{est_ar_yw} and \code{est_ar_dlw} are mainly
#' intended as "internal helper" functions.
#' Therefore, these functions do not check the validity of the input parameters.
#'
#' @param obj either a "time series" object (i.e \code{as.matrix(obj)}
#'     returns an \eqn{(N,m)}-dimensional numeric matrix)
#'     or an \code{\link{autocov}} object which represents an (estimated) autocovariance function.
#'     The type of the \code{autocov} object is irrelevant since \code{est_ar} always uses the
#'     slot \code{obj$gamma} which contains the autocovariance function.
#' @param p.max (integer or \code{NULL}) Maximum order of the candidate AR models.
#'     For the default choice see below.
#' @param p.min (non negative integer) Minimum order of the candidate AR models.
#'     Only used by \code{est_ar_ols}.
#' @param ic (character string) Which information criterion shall be used to find the optimal order.
#'     Note that \code{ic="max"} means that an AR(p) model with \code{p=p.max} is estimated.
#'     Default is \code{ic="AIC"}.
#' @param penalty is a scalar (or NULL) which determines the "penalty" per parameter of the model.
#'     Note that this parameter (if not \code{NULL}) overrides the paramater \code{ic}.
#' @param method Character string giving the method used to fit the model.
#'     Note that 'yule-walker' and 'durbin-levinson-whittle' are (up to numerical errors) equivalent
#'     and that the choice 'ols' is only available for a "time-series" object \code{obj}.
#' @param mean_estimate Character string giving the method used to estimate the mean \eqn{\mu}.
#'     Default is \code{mean_estimate = "sample.mean"}.
#'     See the details below.
#' @param n.obs Optional integer which gives the sample size \eqn{N}.
#'     This parameter is only used, when \code{obj} is an \code{autocov} object.
#'     If \code{n.obs=NULL} then the slot \code{obj$n.obs} is used.
#'     Note that \code{obj$n.obs=NULL} or \code{obj$n.obs=Inf} refers to the case of a population
#'     autocovariance function, i.e. \eqn{N=\infty}.
#'     \cr
#'     For a "time series" object the sample size is of course set to the number of observations,
#'     i.e. \code{n.obs = nrow(as.matrix(obj))}.
#'     \cr
#'     The sample size \eqn{N} controls the computation of the default maximum order \code{p.max} and
#'     the computation of the information criterion.
#' @param gamma \eqn{(m,m,lag.max+1)}-dimensional array, which contains the (sample) autocovariance function.
#' @param y     \eqn{(N,m)}-dimensional matrix, which contains the sample.
#'
#'
#' @return The function \code{est_ar} returns a list with components
#'     \item{model}{\code{\link{armamod}} object which represents the estimated AR model.}
#'     \item{p}{optimal model order.}
#'     \item{stats}{(p.max+1,4) dimensional matrix which stores the \eqn{\ln\det(\Sigma_p)}{ln det (\Sigma[p])}
#'                   values, the number of parameters and the IC values. See the details below.}
#'     \item{y.mean}{estimate of the mean \eqn{\mu}.}
#'     \item{ll}{The log likelihood of the estimated model.}
#' The "helper" functions \code{est_ar_yw}, \code{est_ar_dlw} and \code{est_ar_ols} return a list with components
#'     \item{a}{\code{(m,m,p)}-dimensional array with the estimated AR coefficients \eqn{a_i}{a[i]}.}
#'     \item{sigma}{\code{(m,m)}-dimensional matrix with the estimated noise covariance \eqn{\Sigma}.}
#'     \item{p}{estimate of the AR order.}
#'     \item{stats}{(p.max+1,4) dimensional matrix which stores the \eqn{\ln\det(\Sigma_p)}{ln det (\Sigma[p])}
#'                   values, the number of parameters and the IC values. See the details below.}
#'     \item{y.mean}{(\code{est_ar_ols} only) estimate of the mean \eqn{\mu}.}
#'     \item{residuals}{(\code{est_ar_ols} only) \code{(n.obs,m)} dimensional matrix with the OLS residuals.}
#'     \item{partial}{(\code{est_ar_dlw} only) \code{(m,m,p.max+1)} dimensional array with
#'                    the (estimated) partial autocorrelation coefficients.}
#'
#' @aliases Yule-Walker Durbin-Levinson-Whittle
#'
#' @export
est_ar = function(obj, p.max = NULL, penalty = NULL, ic = c('AIC','BIC','max'),
                  method = c('yule-walker', 'ols', 'durbin-levinson-whittle'),
                  mean_estimate = c('sample.mean', 'intercept','zero'), n.obs = NULL) {
  
  method = match.arg(method)
  ic = match.arg(ic)
  mean_estimate = match.arg(mean_estimate)
  
  if (inherits(obj, 'autocov')) {
    if (method == 'ols') {
      stop('for method="ols" the input parameter "obj" must be a matrix, time series or data-frame')
    }
    gamma = obj$gamma
    m = dim(gamma)[1] # output dimension
    lag.max = dim(gamma)[3] - 1
    if ( (m*(lag.max+1)) == 0 ) stop('"obj" contains no data')
    names = obj$names
    label = obj$label
    y.mean = rep(NA_real_, m)
    
    # check n.obs
    if (is.null(n.obs)) {
      n.obs = obj$n.obs
      if (is.null(n.obs)) n.obs = Inf
    }
    n.obs = as.numeric(n.obs)[1]
    if (is.finite(n.obs)) {
      n.obs = as.integer(n.obs)[1]
      if (n.obs <= 0) stop('the sample size "n.obs" must be non negative')
    } else {
      n.obs = Inf
    }
    
    if (is.null(p.max)) {
      p.max = max(0, min(12, lag.max, 10*log10(n.obs), floor((n.obs-1)/(m+1))))
    }
    p.max = as.integer(p.max)[1]
    if (p.max < 0) stop('p.max must be a non negative integer!')
  } else {
    y = try(as.matrix(obj))
    if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
      stop('input "obj" must be an "autocov" object or a "time series" ',
           'object which may be coerced to a matrix with "as.matrix(y)"')
    }
    m = ncol(y)      # number of "outputs"
    n.obs = nrow(y)  # sample size (optional parameter n.obs is ignored)
    if (m*n.obs == 0) stop('"obj" contains no data')
    
    names = colnames(y)
    label = NULL
    
    if (is.null(p.max)) {
      p.max = max(0, min(12, 10*log10(n.obs), floor((n.obs-1)/(m+1))))
    }
    p.max = as.integer(p.max)[1]
    if (p.max < 0) stop('p.max must be a non negative integer!')
    
    if (method != 'ols') {
      # compute ACF and mean
      demean = (mean_estimate != 'zero')
      gamma = autocov(y, lag.max = p.max, type = 'covariance',
                      na.action = stats::na.fail, demean = demean)$gamma
      if (demean) {
        y.mean = colMeans(y)
      } else {
        y.mean = double(m)
      }
    }
  }
  
  # set penalty
  if (is.null(penalty)) {
    if (ic == 'max') penalty = -1
    if (ic == 'AIC') penalty = 2/n.obs
    if (ic == 'BIC') penalty = log(n.obs)/n.obs
  }
  penalty = as.vector(penalty)[1]
  
  # now call the estimation methods
  if (method == 'ols') {
    out = est_ar_ols(y, p.max = p.max, mean_estimate = mean_estimate, penalty = penalty)
    y.mean = out$y.mean
  }
  if (method == 'yule-walker') {
    out = est_ar_yw(gamma, p.max = p.max, penalty = penalty)
  }
  if (method == 'durbin-levinson-whittle') {
    out = est_ar_dlw(gamma, p.max = p.max, penalty = penalty)
  }
  
  model = armamod(lmfd(a = dbind(d = 3, diag(m), -out$a)),
                  sigma_L = t(chol(out$sigma)), names = names, label = label)
  
  # log Likelihood
  p = out$p
  ll = unname(out$stats[p+1, 'lndetSigma'])
  if (is.finite(n.obs)) {
    ll = (-1 / 2) * (m*log(2*pi) + m + ll)
  } else {
    ll = NA_real_
  }
  
  return(list(model = model, p = p, stats = out$stats, y.mean = y.mean, ll = ll))
}


#' @rdname est_ar
#' @export
est_ar_yw = function(gamma, p.max = (dim(gamma)[3]-1), penalty = -1) {
  # no input check!
  
  m = dim(gamma)[1]   # number of "outputs"
  
  # G is the covariance matrix of (y[t-p.max]',y[t-p+1]',...,y[t]')'
  G = btoeplitz(C = gamma[,,1:(p.max+1),drop=FALSE])
  
  # cholesky decomposition of G
  R = try(chol(G))
  if (inherits(R,'try-error')) stop('Block Toeplitz matrix G = (gamma[i-j]) is not positive definite!')
  
  # vector of log likelihood values and information criteria
  # compute log(det(sigma_p)), where sigma_p is the noise covariance matrix of the AR(p) model
  # note that R[p+1,p+1] is the cholesky decomposition of sigma_p
  ldS = 2*apply(matrix(log(diag(R)), nrow = m, ncol = p.max+1), MARGIN = 2, FUN = sum)
  
  stats = matrix(NA_real_, nrow = p.max+1, 4)
  colnames(stats) = c('p', 'n.par', 'lndetSigma', 'ic')
  stats[, 'p'] = 0:p.max
  stats[, 'lndetSigma'] = ldS
  stats[, 'n.par'] = stats[,'p']*(m^2)
  stats[, 'ic'] = stats[, 'lndetSigma'] + stats[, 'n.par']*penalty
  
  # optimal order
  p = unname(which.min(stats[, 'ic']) - 1)
  
  # compute/estimate AR model of order p
  if (p > 0) {
    # AR coefficients a = (a[p],...,a[1]) are determined by R[1:p,1:p] t(a) = R[1:p,p+1]
    a = backsolve(R[1:(p*m), 1:(p*m), drop=FALSE], R[1:(m*p), (m*p+1):(m*(p+1)), drop = FALSE])
    # we have to reshuffle the coefficients
    a = t(a)
    dim(a) = c(m, m, p)
    a = a[ , , p:1, drop = FALSE]
    
    sigmaR = R[(p*m+1):((p+1)*m), (p*m+1):((p+1)*m), drop=FALSE]
    sigma = t(sigmaR) %*% sigmaR
  } else {
    a = array(0, dim = c(m, m, 0))
    sigma = matrix(gamma[,,1], nrow = m)
  }
  
  return(list(a = a, sigma = sigma, p = p, stats = stats))
}


#' @rdname est_ar
#' @export
est_ar_dlw = function(gamma, p.max = (dim(gamma)[3]-1), penalty = -1) {
  # no input checks!
  
  m = dim(gamma)[1]
  g0 = matrix(gamma[,,1], nrow = m, ncol = m)  # lag zero covariance
  
  # convert gamma to (m*p.max,m) matrix [gamma(1),...,gamma(p.max)]' = E [x[t-1]',...,x[t-p.max]']' x[t]'
  g.past = gamma[,,iseq(2, p.max+1),drop=FALSE]
  dim(g.past) = c(m, m*p.max)
  g.past = t(g.past)
  # print(g.past)
  
  # convert gamma to (m*p.max,m) matrix [gamma(1)',...,gamma(p.max)']' = E [x[t+1]',...,x[t+p.max]']' x[t]'
  # g.future = aperm(gamma[,,iseq(2,p.max+1),drop=FALSE],c(1,3,2))
  # dim(g.future) = c(m*p.max, m)
  # print(g.future)
  
  # initialize ( <=> order p=0 model)
  a = matrix(0, nrow = m, ncol = m*p.max) # matrix with the coefficients of the forward model
  a0 = a
  b = a            # matrix with the coefficients of the backward model
  sigma = g0       # covariance of the forecasting errors
  sigma0 = sigma
  omega = sigma    # covariance of the backcasting errors
  omegae0 = omega
  
  # partial autocorrelation coefficients
  partial = array(0, dim = c(m, m, p.max+1))
  su = matrix(sqrt(diag(sigma)), nrow = m, ncol = m)
  partial[,,1] = sigma / (su * t(su))
  
  stats = matrix(NA_real_, nrow = p.max+1, 4)
  colnames(stats) = c('p', 'n.par', 'lndetSigma', 'ic')
  stats[, 'p'] = 0:p.max
  stats[, 'n.par'] = stats[,'p']*(m^2)
  
  # optimal model so far
  stats[1, 'lndetSigma'] = log(det(sigma))
  stats[1, 'ic'] = stats[1, 'lndetSigma'] + stats[1, 'n.par']*penalty
  ic.opt = stats[1, 'ic']
  p.opt = 0
  a.opt = a
  sigma.opt = sigma
  
  # i is a matrix of indices, such that we can easily access the contents of g.future, g.past, a, b, ...
  # e.g. i[,1] = 1:m, i[,2] = (m+1):2*m, ...
  i = matrix(iseq(1,m*p.max), nrow = m, ncol = p.max)
  
  for (p in iseq(1, p.max)) {
    # compute order p model
    
    sigma0 = sigma
    omega0 = omega
    if (p > 1) a0[, i[,1:(p-1)]] = a[, i[,1:(p-1)]]
    
    # update forward model
    # E v[t-p] y[t]', where v[t-p] = y[t-p] - b[1] y[t-p+1] - ... - b[p-1] y[t-1]
    # is the backcasting error
    if (p > 1) {
      Evy = g.past[i[, p], , drop=FALSE] - b[,i[,1:(p-1)],drop=FALSE] %*% g.past[i[,(p-1):1],,drop=FALSE]
    } else {
      Evy = g.past[i[, p], , drop=FALSE]
    }
    aa = t(solve(omega, Evy))
    sigma = sigma - aa %*% Evy
    a[, i[, p]] = aa
    if (p>1) {
      a[, i[, 1:(p-1)]] = a[, i[, 1:(p-1)], drop=FALSE] - aa %*% b[, i[, (p-1):1], drop = FALSE]
    }
    
    # update backward model
    # E u[t] y[t-p]' = t( E v[t-p] y[t]' ), where u[t] = y[t] - a[1] y[t-1] - ... - a[p-1] y[t-p+1]
    # is the forecasting error
    # Euy = g.future[i[,p],,drop=FALSE] - a0[,i[,iseq(1,p-1)],drop=FALSE] %*%
    #       g.future[i[,rev(iseq(1,p-1))],,drop=FALSE]
    Euy = t(Evy)
    bb = t(solve(sigma0, Euy))
    omega = omega - bb %*% Euy
    b[, i[, p]] = bb
    if (p>1) {
      b[, i[, 1:(p-1)]] = b[, i[, 1:(p-1)], drop=FALSE] - bb %*% a0[, i[, (p-1):1], drop=FALSE]
    }
    
    # partial autocorrelations
    su = matrix(sqrt(diag(sigma0)), nrow=m, ncol=m)
    sv = matrix(sqrt(diag(omega0)), nrow=m, ncol=m,byrow = TRUE)
    # print(cbind(Euy,sigma0,omega0,su,sv))
    partial[,,p+1] = Euy / (su * sv)
    
    stats[p+1, 'lndetSigma'] = log(det(sigma))
    stats[p+1, 'ic'] = stats[p+1, 'lndetSigma'] + stats[p+1, 'n.par']*penalty
    
    # cat(p,'\n')
    # print(Evy)
    # print(omega)
    # print(eigen(omega)$values)
    # print(sigma)
    # print(eigen(sigma)$values)
    # print(ldS)
    # print(ic)
    if (stats[p+1, 'ic'] < ic.opt) {
      # update optimal model
      ic.opt = stats[p+1, 'ic']
      a.opt = a
      sigma.opt = sigma
      p.opt = p
    }
  }
  
  a = a.opt[, iseq(1,m*p.opt), drop = FALSE]
  dim(a) = c(m, m, p.opt)
  return(list(a = a, sigma = sigma.opt, p = p.opt, stats = stats, partial = partial))
}

#' @rdname est_ar
#' @export
est_ar_ols = function(y, p.max = NULL, penalty = -1,
                      mean_estimate = c('sample.mean', 'intercept','zero'), p.min = 0L) {
  # only some basic input checks
  mean_estimate = match.arg(mean_estimate)
  intercept = (mean_estimate == 'intercept')
  
  # coerce data objects to matrix
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  m = ncol(y)      # number of "outputs"
  n.obs = nrow(y)  # sample size
  if ( (m*n.obs == 0) ) stop('"y" contains no data')
  
  p.min = as.integer(p.min)[1]
  if (p.min < 0) stop('minimum order "p.min" must be a non negative integer')
  
  if (is.null(p.max)) {
    p.max = max(p.min, min(12, 10*log10(n.obs)/m, floor((n.obs-1)/(m+1)) ))
  }
  p.max = as.integer(p.max)[1]
  if ( (p.max < 0) ) stop('maximum order "p.max" must be a non negative integer!')
  if ( (p.max < p.min) ) stop('maximum order "p.max" is smaller than the minimum order "p.min"!')
  if ( (n.obs-p.max) < (m*p.max + intercept) ) {
    stop('sample size N is not sufficient for the desired maximum order "p.max"')
  }
  
  y.mean = double(m)
  if (mean_estimate == 'sample.mean') {
    y.mean = colMeans(y) # sample mean
    y = y - matrix(y.mean, nrow = n.obs, ncol = m, byrow = TRUE)
  }
  
  # create a "big" data matrix X, where the rows are of the form
  #    (y[t,], y[t-1,],...,y[t-p.max,])
  
  # use "btoeplitz"
  # dim(y) = c(n.obs,n,1)
  # y = aperm(y, c(3,2,1))
  # junk = array(0,dim=c(1,n,p.max+1))
  # junk[,,1] = y[,,1]
  # X = btoeplitz(R = junk, C = y)
  
  # use "for loop"
  X = matrix(NA_real_, nrow = n.obs, ncol = m*(p.max+1))
  for (i in (0:p.max)) {
    X[(1+i):n.obs, (m*i+1):(m*(i+1))] = y[1:(n.obs-i),]
  }
  # print(X[1:(n*(p.max+2)),])
  
  stats = matrix(NA_real_, nrow = p.max-p.min+1, 4)
  colnames(stats) = c('p', 'n.par', 'lndetSigma', 'ic')
  stats[, 'p'] = p.min:p.max
  ic.opt = Inf
  u.opt = matrix(NA_real_, nrow = n.obs, ncol = m)
  p.opt = NA_integer_
  
  # estimate AR models with order p = p.min, 1, ..., p.max
  for (p in (p.min:p.max)) {
    if (p == 0) {
      # AR(0) model
      a = matrix(0, nrow = m, ncol = 0)
      if (intercept) {
        d = colMeans(y)
        u = y - matrix(d, nrow = n.obs, ncol = m, byrow = TRUE)
      } else {
        d = double(m)
        u = y
      }
    } else {
      # estimate coefficients by OLS
      out = stats::lsfit(X[(1+p):n.obs, (m+1):(m*(p+1)), drop=FALSE],
                         X[(1+p):n.obs, 1:m, drop=FALSE], intercept = intercept)
      a = t(out$coef)
      if (intercept) {
        d = a[, 1]
        a = a[, -1, drop = FALSE]
      } else {
        d = NULL
      }
      u = out$residuals
    }
    sigma = crossprod(u)/(n.obs - p)
    # log(det(sigma)) may generate NaN's, suppress warning
    i = which(stats[,'p'] == p)
    stats[i, 'lndetSigma'] = suppressWarnings( log(det(sigma)) )
    stats[i, 'n.par'] = (m^2)*p + intercept*m
    stats[i, 'ic'] = stats[i, 'lndetSigma'] + stats[i, 'n.par']*penalty
    # take care of NaN's
    if ( is.finite(stats[i, 'ic']) && (stats[i, 'ic'] < ic.opt)) {
      # we have found a better model!
      a.opt = a
      d.opt = d
      sigma.opt = sigma
      p.opt = p
      ic.opt = stats[i, 'ic']
      u.opt[1:p, ] = NA_real_
      u.opt[(p+1):n.obs, ] = u
    }
  }
  
  # this should not happen
  if (is.na(p.opt)) stop('could not find an optimalt AR model')
  
  # coerce a.opt to 3-D array
  dim(a.opt) = c(m, m, p.opt)
  
  if (intercept) {
    # estimated mean from intercept mu = (I -a[1] - ... - a[p])^{-1} d
    if (p.opt == 0) {
      y.mean = d.opt
    } else {
      # a(1) = I - a[1] - ... - a[p]
      # print( apply(a.opt, MARGIN = c(1,2), FUN = sum) )
      a1 = diag(m) - matrix(apply(a.opt, MARGIN = c(1,2), FUN = sum), nrow = m, ncol = m)
      y.mean = try(solve(a1, d.opt))
      if (inherits(a1, 'try-error')) {
        # model has a unit root a(1) is singular!!!!
        y.mean = rep(NaN, m)
      }
    }
  }
  
  # # log likelihood
  # ll = unname((-(n.obs - p.opt)/2)*(m*log(2*pi) + m + stats[p.opt+1,'lndetSigma']))
  # # cat('est_ar_ols', (n.obs-p.opt)/2, m*log(2*pi), m, stats[p.opt+1, 'lndetSigma'], ll, '\n')
  
  # return the best model!
  return(list(a = a.opt, sigma = sigma.opt, p = p.opt,
              stats = stats, y.mean = y.mean, residuals = u.opt))
}




