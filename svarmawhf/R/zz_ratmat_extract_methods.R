# subsetting for rational matrices ###############################################  

# helper function:
# check the result of subsetting x[]
# based on the result of subsetting an "ordinary" matrix 
extract_matrix_ = function(m, n, n_args, is_missing, i, j) {
  idx = matrix(iseq(1, m*n), nrow = m, ncol = n)
  
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(idx)
  
  if (n_args == 1) {
    # x[i]
    idx = idx[i]
    if (any(is.na(idx))) stop('index out of bounds')
    idx = matrix(idx, nrow = length(idx), ncol = 1)
    return(idx)
  }
  
  if (is_missing[1]) {
    # x[,j]
    return(idx[, j, drop = FALSE])
  }
  if (is_missing[2]) {
    # x[i,]
    return(idx[i, , drop = FALSE])
  }
  # x[i,j]
  return(idx[i, j, drop = FALSE])
}

#' Extract Parts of a Rational Matrix
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' The subsetting operation \code{x[,]} for rational matrices works analogously to the subsetting of 
#' ordinary matrices. However, this operator is not implemented for \code{\link{lmfd}} objects. 
#' (x[i,j] throws an error if x is an \code{lmfd} object!) The result is an object of the same class. 
#' \cr
#' The \code{$} operator may be used to extract the \eqn{a}, \eqn{b} polynomial of a left matrix 
#' fraction description (\code{\link{lmfd}} object) and the matrices 
#' \eqn{A,B,C,D} of a state space representation (\code{\link{stsp}} object).
#' 
#' \itemize{
#' \item \code{x[]} or \code{x[,]} simply return the original object.
#' \item \code{x[i]} returns a "vector", i.e. an \eqn{(s,1)} dimensional matrix. 
#' \item \code{x[i,]},\code{x[,j]} or \code{x[i,j]} return a rational matrix with rows selected by \code{i} 
#'       and columns selected by \code{j}. 
#' }
#' Note that "named" arguments are not supported (in order to simplify the coding). 
#' 
#' In order to have a finer control, one may e.g. use \code{unclass(x)[,,]}.
#' 
#'
#' @param x a rational matrix, i.e. a \code{\link{polm}}, \code{\link{lmfd}}, \code{\link{stsp}}, 
#'          \code{\link{pseries}} or \code{\link{zvalues}} object.
#' @param i,j indices (integer or boolean vector) 
#' @param name character (A,B,C,D for statespace objects and a,b for LMFD objects)
#'
#' @return rational matrix
#' 
#' @rdname extract
#' @name extract
#' @export
#' 
#' @examples 
#' # polynomial matrices 
#' a = test_polm(dim = c(3,2), degree = 1)
#' a[]          # returns a
#' a[,]         # returns a
#' a[c(1,3,6)]  # returns a "vector" with the (1,1), (3,1) and (3,2) element of a
#' a[1,]        # returns the first row of a 
#' a[,2]        # returns the second column of a 
#' a[c(TRUE,FALSE,TRUE),c(FALSE, TRUE)] # returns a 2 by 1 matrix 
#' a[c(1,1),c(2,1)] # returns a 2 by 2 matrix 
#' 
#' \dontrun{
#' # the subsetting operator [,] is not implemented for "lmfd" objects
#' test_lmfd(dim = c(2,2), degrees = c(1,1))[1,1]
#' 
#' a[i=1, j=2] # throws an error, since "named" arguments are not allowed.
#' }
'[.polm' = function(x,i,j) {
  # print(sys.call())
  if (!is.null(names(sys.call()))) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 1
  is_missing = c(missing(i), missing(j))
  # x[] or x[,]
  if ( ((n_args==1) && is_missing[1]) || all(is_missing) ) return(x)
  
  x = unclass(x)
  d = dim(x)
  m = d[1]
  n = d[2]
  p = d[3] - 1
  
  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)
  
  if (length(idx) == 0) {
    # result is an "empty" rational matrix
    x = polm(array(0, dim = c(nrow(idx), ncol(idx), 1)))
    return(x)
  }
  
  if (n_args == 1) { 
    dim(x) = c(m*n,1,p+1)
    x = polm(x[i,1,,drop = FALSE])
    return(x)
  }
  
  if (is_missing[1]) {
    # x[,j]
    x = polm(x[,j,,drop = FALSE])
    return(x)
  }
  
  if (is_missing[2]) {
    # x[i,]
    x = polm(x[i,,,drop=FALSE])
    return(x)
  }
  
  # x[i,j]
  x = polm(x[i,j,,drop=FALSE])
  return(x)
}

#' @rdname extract
#' @export
'[.lmfd' = function(x,i,j) {
  stop('The subsetting operator [] is not implemented for "lmfd" objects.')
}
  


#' @rdname extract
#' @export
'$.lmfd' = function(x, name) {
  i = match(name, c('a','b'))
  if (is.na(i)) stop('reference to "',name, '" is not defined')
  d = attr(x,'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  x = unclass(x)
  if (i == 1) {
    return(polm(array(x[,iseq(1,m*(p+1))], dim = c(m,m,p+1))))
  }
  if (i == 2) {
    return(polm(array(x[,iseq(m*(p+1)+1,m*(p+1)+n*(q+1))], dim = c(m,n,q+1))))
  }
  # this should not happen
  stop('unknown reference')
}

#' Replace Parts of a Polynomial Matrix
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' The assigment operation \code{x[,] <- value} for polynomial matrices works quite analogously 
#' to the assigment operation of "ordinary" matrices. 
#' 
#' Note that "named" arguments are not supported (in order to simplify the coding). 
#'
#' @param x \code{\link{polm}} object
#' @param i,j indices
#' @param value  Either a \code{\link{polm}} object, or a vector/matrix/array which may be coerced to 
#'               a \code{polm} object by \code{polm(value)}.
#'
#' @export
#' @rdname replace
#' @name replace
#'
#' @examples
#' a = test_polm(dim = c(3,2), degree = 1)
#' a
#' a[FALSE] = 0   # no items to replace, a is not changed
#' a
#' a[lower.tri(matrix(0, nrow = 3, ncol = 2))] = 0 # set elements below the diagonal equal to zero
#' a
#' a[3,1] = c(1,-1) # set (3,1) element
#' a
#' a[1:2, 2:1] = c(0,1)
#' a
#' a[2, ] = test_polm(dim = c(1,2), degree = 4)
#' a
#' a[, 1] = test_polm(dim = c(2,1), degree = 4) # this gives a warning
#' a
#' 
#' \dontrun{
#' a[i=1] = 0   # named arguments are not supported
#' }  
"[<-.polm" = function(x,i,j,value) {
  names_args = names(sys.call())
  # print(sys.call())
  # print(names_args)
  if (!all(names_args[-length(names_args)] == '')) {
    stop('named dimensions are not supported')
  }
  n_args = nargs() - 2
  # print(n_args)
  is_missing = c(missing(i), missing(j))
  # print(is_missing)
  
  x = unclass(x)
  dx = dim(x)
  m = dx[1]
  n = dx[2]
  p = dx[3] - 1 
  
  idx = try(extract_matrix_(m, n, n_args, is_missing, i, j), silent = TRUE)
  if (inherits(idx, 'try-error')) stop('index/subscripts out of bounds')
  # print(idx)
  idx = as.vector(idx)
  # print(length(idx))

  if (!inherits(value,'polm')) {
    # coerce right hand side 'value' to polm object  
    value = try( polm(value) )
    if (inherits(value, 'try-error')) {
      stop('Could not coerce the right hand side to a polm object!')
    }
  }
  value = unclass(value)
  dv = dim(value)
  
  # no items to replace: return original object
  if (length(idx) == 0) return(polm(x))
  
  if ((dv[1]*dv[2]) == 0) stop('replacement has length 0')
  
  # bring degrees of 'x' and of 'value' in line
  if (dv[3] > dx[3]) {
    x = dbind(d = 3, x, array(0, dim = c(dx[1], dx[2], dv[3] - dx[3])) )
    p = dv[3] - 1
  }
  if (dv[3] < dx[3]) {
    value = dbind(d = 3, value, array(0, dim = c(dv[1], dv[2], dx[3] - dv[3])) )
  }
  
  # coerce 'x' and 'value' to "vectors"
  dim(x) = c(m*n, p+1)
  dim(value) = c(dv[1]*dv[2], p+1)
  # print(value)
  
  # extend 'value' if needed
  if ( (length(idx) %% nrow(value)) != 0 ) {
    warning('number of items to replace is not a multiple of replacement length')
  }
  value = value[rep_len(1:nrow(value), length(idx)),,drop = FALSE]
  # print(value)
  
  # plug in new values
  x[idx,] = value
  # reshape
  dim(x) = c(m, n, p+1)
  # re-coerce to polm
  x = polm(x)
  
  return(x) 
}
