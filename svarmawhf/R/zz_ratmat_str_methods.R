# str.____ methods ##############################################################

#' Display the Structure of Objects
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' 
#' @param object an object
#' @param ... not used
#'
#' @return invisible(NULL)
#' 
#' @rdname str
#' @name str methods
#' @export
str.polm = function(object, ...) {
  d = dim(unclass(object))
  cat('polynomial matrix [',d[1],',',d[2],'] with degree <= ',d[3]-1,'\n', sep = '')
  return(invisible(NULL))  
}

#' @rdname str
#' @export
str.lmfd = function(object, ...) {
  d = attr(object, 'order')
  cat('left matrix fraction description [',d[1],',',d[2],'] with degrees (p = ', 
      d[3], ', q = ', d[4],')\n', sep = '')
  return(invisible(NULL))  
}
