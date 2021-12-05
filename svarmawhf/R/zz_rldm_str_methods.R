# str.____ methods ##############################################################

#' Display the Structure of Objects
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' 
#' @param object an object
#' @param ... not used
#'
#' @return invisible(NULL)
#'
#' @rdname str
#' @name str methods
#' @export
str.armamod = function(object, ...) {
  d = attr(object$sys, 'order')
  cat('ARMA model [',d[1],',',d[2],'] with orders p = ', d[3], ' and q = ', d[4], '\n', sep = '')
  return(invisible(NULL))
}
