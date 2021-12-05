# dim.____ methods ##############################################################

#' Dimensions of Objects
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' 
#' @param x object
#'
#' @return returns a named integer vector
#' \item{\code{c(m,n,p)}}{for polynomial matrices (i.e. for \code{\link{polm}} objects).}
#' \item{\code{c(m,n,p,q)}}{for left matrix fraction descriptions (i.e. for \code{\link{lmfd}} objects).}
#' \item{\code{c(m,n,s)}}{for state space representations (i.e. for \code{\link{stsp}} objects).}
#' \item{\code{c(m,n,lag.max)}}{for impulse response functions (i.e. for \code{\link{pseries}} objects).}
#' \item{\code{c(m,n,n.f)}}{for frequency response functions (i.e. for \code{\link{zvalues}} objects).}
#' 
#' @rdname dim
#' @name dim methods
#' @export
dim.polm = function(x) {
  d = dim(unclass(x))
  d[3] = d[3] - 1
  names(d) = c('m','n','p')
  return(d)
}

#' @rdname dim
#' @export
dim.lmfd = function(x) {
  d = attr(x, 'order')
  names(d) = c('m','n','p','q')
  return(d)
}


#' @rdname dim
#' @export
dim.pseries = function(x) {
  d = dim(unclass(x))
  d[3] = d[3] - 1
  names(d) = c('m','n','lag.max')
  return(d)
}