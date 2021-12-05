# Misc Tools and Utilities ###############################################

#' Sequence generation
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#'
#' As opposed to the standard \code{\link{seq}} command the function \code{iseq} returns an empty
#' vector if the starting value is larger than the end value.
#'
#' More general than \link{seq_len} because the sequence does not need to start from 1.
#'
#' @param from,to the starting and end values of the sequence.
#'
#' @return \code{iseq} returns an empty integer vector if \code{to} is less than \code{from} and
#'    \code{seq(from,to)} else.
#' @export
#'
#' @examples
#' iseq(0,1) # => c(0,1)
#' iseq(1,0) # => integer(0)
#' seq(1,0)  # => c(1,0)
iseq = function(from = 1, to = 1) {
  if (to<from) return(integer(0))
  return( seq(from = from, to = to) )
}


# expand_letters is an internal function
# #' Creating letters
# #'
# #' Create a 'character' vector of length n out of the letters 'l' if 'l' has less than n entries,
# #' then combinations are created.
# #'
# #' @note
# #' Should eventually be internal.
# #'
# #' @param n Integer. Number of letters to be created.
# #' @param l Vector of characters. Must be at least of length 2.
# #'
# #' @return Vector of characters of length n.
# #' @export
# #'
# #' @examples
# #' expand_letters(4, l = c('a','b'))
# #' expand_letters(6, l = c('a','b'))
expand_letters = function(n, l = letters) {
  if (n == 0) return(character(0))
  if ((n>1) && (length(l) <= 1)) stop('"l" must have at least 2 entries!')
  l0 = l
  while (length(l) < n) {
    l = outer(l,l0,FUN = function(a,b) {paste(b,a,sep='')})
  }
  l[1:n]
}


#' Create Test Array
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' The helper function \code{test_array} creates arrays of given dimension. 
#' 
#' The function \code{test_array} creates an array with either 
#' \itemize{
#' \item entries of the form \code{x[i,j,...] = i*10^(n-1) + j*10^(n-2) + ...},
#'       where \code{n} is the number of dimensions, 
#' \item or randomly generated (standard normal) entries.
#' }
#' 
#' This function is mainly used to test operations on arrays
#' (like \code{\link{btoeplitz}}, \code{\link{bhankel}}, \code{\link{bmatrix}} and 
#' \code{\link{dbind}}).
#'
#' @param dim Integer vector.
#' @param random Boolean. Default set to FALSE.
#'   If TRUE, the elements are randomly generated (standard normal).
#' @param dimnames (logical) decides whether \code{dimnames} should be created. 
#'
#' @return Real valued array of dimension \code{dim}.
#' @export
#'
#' @examples
#' test_array(dim = 5)
#' test_array(dim = c(2,4), dimnames = TRUE)
#' test_array(dim = c(2,4,3), dimnames = FALSE)
#'
#' \dontrun{
#' # the examples below throw an error
#' test_array(dim = c())
#' test_array(dim = c(2,-1,2))
#' }
test_array = function(dim, random = FALSE, dimnames = FALSE) {
  # check input parameter "dim"
  dim = as.vector(as.integer(dim))
  if ((length(dim)==0) || (min(dim) < 0)) {
    stop('"dim" must be a (non empty) vector of non negative integers!')
  }
  n.dims = length(dim)
  
  if (min(dim)==0) {
    # empty array
    x = array(0, dim=dim)
  } else {
    if (random) {
      x = array(stats::rnorm(prod(dim)), dim = dim)
    } else {
      x = (1:dim[1])
      for (i in (iseq(2,n.dims))) {
        x = outer(10*x,(1:dim[i]),FUN = '+')
      }
      x = array(x, dim = dim)
    }
  }
  
  # set dimnames
  if (dimnames) {
    dimnames.x = as.list(1:n.dims)
    names(dimnames.x) = expand_letters(n.dims, l = LETTERS)
    for (i in (1:n.dims)) {
      if (dim[i]>0) {
        dimnames.x[[i]] = paste(names(dimnames.x)[i],1:dim[i],sep = '=')
      } else {
        dimnames.x[[i]] = character(0)
      }
    }
    dimnames(x) = dimnames.x
  }
  
  return(x)
}



# Array Tools ##################################################################

#' Bind Arrays
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' \code{dbind(d, x, y, ...)} concatenates/binds an arbitrary number of arrays \code{x, y,...} along
#' the dimension \code{d}.
#' For matrices, \code{dbind(d = 1, x, y, ...)} is (essentially) equivalent to
#' \code{\link{rbind}} and \code{dbind(d = 2, x, y, ...)} corresponds to \code{\link{cbind}}.
#' If the number of dimensions of an argument, \code{x} say,  is less than \code{d}, then this argument
#' is treated as an array of dimension \code{c(dim(x),1,..,1)}.
#'
#' The procedure makes some effort to keep the \code{dimnames} attribute of the arguments.
#'
#' @param d Integer. Concatenate arrays along the dimension \code{d}.
#' @param ... arrays.
#'
#' @return Array
#' @export
#'
#' @examples
#' x = test_array(dim = c(2,3,1))
#' y = test_array(dim = c(2,3,1), dimnames = TRUE)
#' z = test_array(dim = c(2,3,3), dimnames = TRUE)
#'
#' # Bind along dimension 1 (row-binding for matrices)
#' dbind(1, x)
#' dbind(1, x, y)
#'
#' # Bind along dimension 2 (col-binding for matrices)
#' dbind(2, x, y)
#'
#' # Bind along dimension 3
#' dbind(3, x, y)
#' dbind(3, x, y, z)
#'
#' # Example that throws an error
#' \dontrun{
#' dbind(1, x, y, z) # throws an error, since the array x,y,z are not compatible
#' }
dbind = function(d = 1, ...) {
  d = as.integer(as.vector(d))[1]
  if (d < 1) stop('"d" must be a positive integer!')
  
  args = list(...)
  n.args = length(args)
  if (n.args == 0) {
    stop('no arrays to bind!')
  }
  
  a = args[[1]]
  if (is.vector(a)) {
    dimnames.a = names(a)
    dim(a) = length(a)
    dimnames(a) = list(dimnames.a)
  }
  dim.a = dim(a)
  dimnames.a = dimnames(a)
  if (d > length(dim.a)) {
    dim(a) = c(dim(a),rep(1, d - length(dim.a)))
    dimnames(a) = c(dimnames.a, vector('list', d - length(dim.a)))
  }
  dim.a = dim(a)
  dimnames.a = dimnames(a)
  
  if (n.args == 1) {
    return(a)
  }
  
  b = args[[2]]
  if (is.vector(b)) {
    dimnames.b = names(b)
    dim(b) = length(b)
    dimnames(b) = list(dimnames.b)
  }
  dim.b = dim(b)
  dimnames.b = dimnames(b)
  if (d > length(dim.b)) {
    dim(b) = c(dim(b),rep(1, d - length(dim.b)))
    dimnames(b) = c(dimnames.b, vector('list', d - length(dim.b)))
  }
  dim.b = dim(b)
  dimnames.b = dimnames(b)
  
  n.dim = length(dim.a)
  if ( (length(dim.b) != n.dim) || ( (n.dim > 1) && (any(dim.a[-d]!=dim.b[-d]))) ) {
    stop('arrays are not compatible!')
  }
  
  # bind arrays a and b -> c
  dim.c = dim.a # dimension of the result
  dim.c[d] = dim.a[d] + dim.b[d]
  p = c((1:n.dim)[-d],d)
  pp = (1:n.dim)
  pp[p] = 1:n.dim
  a = aperm(a, p)
  b = aperm(b, p)
  dim.c = dim.c[p]
  c = array(c(a,b), dim = dim.c)
  c = aperm(c, pp)
  
  # take care of dimnames
  if ( (is.null(dimnames.a)) && (!is.null(dimnames.b)) ) {
    dimnames.a = dimnames.b
    dimnames.a[[d]] = character(0)
  }
  if ( (is.null(dimnames.b)) && (!is.null(dimnames.a)) ) {
    dimnames.b = dimnames.a
    dimnames.b[[d]] = character(0)
  }
  
  if ( !is.null(dimnames.a) ) {
    dimnames.c = dimnames.a
    if (is.null(names(dimnames.c))) names(dimnames.c) = names(dimnames.b)
    
    for (i in (1:n.dim)) {
      if ( (names(dimnames.c)[i] == '') && (!is.null(names(dimnames.b))) ) {
        names(dimnames.c)[i] = names(dimnames.b)[i]
      }
      
      if (i != d) {
        if ( (length(dimnames.c[[i]]) == 0) && (length(dimnames.b[[i]]) > 0) ) {
          dimnames.c[[i]] = dimnames.b[[i]]
        }
      } else {
        if ( (length(dimnames.a[[i]]) == 0) || (length(dimnames.b[[i]]) == 0) ) {
          dimnames.c[[i]] = character(0)
        } else {
          dimnames.c[[i]] = c(dimnames.a[[i]], dimnames.b[[i]])
        }
      }
      
    }
    dimnames(c) = dimnames.c
  }
  
  if (n.args == 2) {
    return(c)
  }
  
  args[[2]] = c
  cc = do.call(dbind, args = c(list(d = d), args[-1]))
  return(cc)
}




#' Coerce arrays to data frames
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' This helper function creates a \code{data.frame} from a given array.
#'
#' @param x Matrix or array.
#' @param rows,cols integer vectors. These two vectors define a partition of the "dimensions" (1,...,n),
#'        where n is the number of dimensions of x (i.e. length(dim(x))).
#'        If either of the two is missing, then the complement is used.
#'        At least one of the arguments "rows" and "cols" has to be given.
#'
#' @return data.frame
#' @export
#'
#' @seealso The helper function \code{array2data.frame} is used
#'          internally for the \code{\link{plot methods}}. 
#'          The function \code{\link{bmatrix}} (which coerces arrays to matrices) is a simplified
#'          version of \code{array2data.frame}.
#'
#' @examples
#' # test array
#' x = test_array(dim = c(2,3,2), dimnames = TRUE)
#' array2data.frame(x, cols = c(1,3,2))
#' array2data.frame(x, rows = 1)
#' array2data.frame(x, rows = 2:1)
#' array2data.frame(x, rows = c(2,1,3))
#'
#' # consider a pseudo socio economic data set
#' x = test_array(dim = c(2,4,5), random = TRUE)
#' dimnames(x) = list(sex=c('female','male'),
#'                    education = c('none','primary','high','university'),
#'                    age = c('<20','30-40','40-50','50-60','>60'))
#' array2data.frame(x, cols = 1)
#' array2data.frame(x, cols = c(1,2))
#'
array2data.frame = function(x, rows = NULL, cols = NULL) {
  
  # check parameter "x"
  if ( (!is.array(x)) || (min(dim(x)) <= 0) || (length(dim(x)) == 0) ) {
    stop('"x" must be a (non empty) array!')
  }
  dim.x = dim(x)
  n.dims = length(dim.x)
  dims = (1:n.dims)
  dimnames.x = dimnames(x)
  if (is.null(dimnames.x)) stop('"x" must have a complete "dimnames" attribute!')
  names.dimnames.x = names(dimnames.x)
  if ( (is.null(names.dimnames.x)) || (any(names.dimnames.x =='')) ||
       (length(unique(names.dimnames.x)) < n.dims ) ) {
    stop('"x" must have a complete "dimnames" attribute!')
  }
  for (i in dims) {
    if ( (is.null(dimnames.x[[i]])) || (any(dimnames.x[[i]] =='')) ||
         (length(unique(dimnames.x[[i]])) < dim.x[i]) ) {
      stop('"x" must have a complete "dimnames" attribute!')
    }
  }
  
  # check input parameters "rows" and "cols"
  if (is.null(rows) && is.null(cols)) {
    stop('missing parameters "rows" and "cols"!')
  }
  
  # check "rows"
  if (!is.null(rows)) {
    rows = as.integer(as.vector(rows))
    # rows must be subset of dim.x
    if ( (length(rows) > 0) && ( (min(rows) < 1) || (max(rows) > n.dims) ) ) {
      stop('parameter "rows" does not correspond to a selection of dimensions (of "x")!')
    }
  }
  
  # check "cols"
  if (!is.null(cols)) {
    cols = as.integer(as.vector(cols))
    # cols must be subset of dim.x
    if ( (length(cols) > 0) && ( (min(cols) < 1) || (max(cols) > n.dims) ) ) {
      stop('parameter "cols" does not correspond to a selection of dimensions (of "x")!')
    }
  }
  
  if (is.null(rows)) {
    if (length(cols) > 0) {
      rows = dims[-cols]
    } else {
      rows = dims
    }
  }
  if (is.null(cols)) {
    if (length(rows) > 0) {
      cols = dims[-rows]
    } else {
      cols = dims
    }
  }
  
  # cat('rows:',rows,'cols:',cols,'\n')
  
  # check that the union of "rows" and "cols" = "dim.x"
  if (!isTRUE(all.equal(dims, sort(c(rows,cols))))) {
    stop('the parameters "rows" and "cols" do not correspond to a partition of the set of dimensions (of "x")!')
  }
  
  # coerce x to a matrix
  x = aperm(x, c(rows,cols))
  dim(x) = c(prod(dim.x[rows]),prod(dim.x[cols]))
  
  if (length(rows) > 0) {
    cases = expand.grid(dimnames.x[rows], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
    # print(cases)
    # print(x)
    x = cbind(cases, x)
    # print(x)
    # print(str(x))
  } else {
    x = data.frame(x)
  }
  
  if (length(cols) > 0) {
    variables = expand.grid(dimnames.x[cols], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    # print(variables)
    # print(str(variables))
    variables = apply(variables, MARGIN = 1, FUN = paste, collapse='.')
    # print(variables)
  } else {
    variables = 'value'
  }
  
  if (length(rows)>0) {
    colnames(x) = c(names.dimnames.x[rows],variables)
  } else {
    colnames(x) = variables
  }
  
  return(x)
}


# Block Matrices #################################################################

#' Block Diagonal Matrix
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Combine two or more matrices to a block-diagonal matrix. The functions supports boolean, integer,
#' numeric and complex matrices (and vectors). The procedure makes some effort to retain the
#' (col- / row-) names of the inputs.
#'
#' @param ... matrices or vectors. Vectors are treated as diagonal matrices.
#'            If no input arguments are provided then \code{bdiag} returns \code{NULL}.
#'
#' @return block diagonal matrix (or \code{NULL} if \code{bdiag} has been called without inputs).
#' @export
#'
#' @examples
#' A = matrix(TRUE, nrow = 2, ncol = 3)
#' colnames(A) = paste('A',1:3,sep ='.')
#' B = rep(2L,3)
#' names(B) = paste('B',1:3,sep = '.')
#' C = matrix(0, nrow = 3, ncol = 0)
#' rownames(C) = paste('C',1:3)
#' D = matrix(4, nrow = 2, ncol = 3)
#' E = matrix(complex(real=5), nrow = 2, ncol = 2)
#' rownames(E) = paste('E',1:2,sep='.')
#' colnames(E) = paste('E',1:2,sep='.')
#'
#' bdiag()
#' bdiag(NULL,NULL)
#'
#' X = bdiag(A,NULL) # NULL arguments are skipped
#' X
#' str(X)            # output is of type 'logi'
#'
#' X = bdiag(A,NULL,B)
#' X
#' str(X)            # output is of type 'int'
#'
#' # note the action of the "empty" (3 times 0) matrix C
#' X = bdiag(A,C,B)
#' X
#' str(X)            # output is of type 'num'
#'
#' X = bdiag(A,B,C,D,E)
#' X
#' str(X)            # output is of type 'cplx'
#'
#' \dontrun{
#' # the inputs must be vectors or matrices. arrays are not supported.
#' bdiag(A, array(1, dim = c(2,3,1)))
#'
#' # character matrices are problematic, since it is not clear how to set the
#' # non diagonal elements. Therefore, the following statement throws an error.
#' bdiag(A, matrix('B', nrow = 2, ncol = 1), 3:4)
#' }
bdiag = function(...) {
  
  args = list(...)
  n_args = length(args)
  
  # no inputs -> return NULL
  if (n_args == 0) return(NULL)
  
  # skip 'NULL' arguments
  i = sapply(args, FUN = function(x) (!is.null(x)))
  
  if (sum(i) == 0) return(NULL)
  args = args[i]
  n_args = length(args)
  
  # check type of arguments
  i = sapply(args, FUN = function(x) ( !( is.logical(x) || is.integer(x) || is.numeric(x) || is.complex(x) ) ) )
  if (any(i)) stop('inputs must be of type logical, integer, numeric or complex!')
  i = sapply(args, FUN = function(x) ( !( is.vector(x) || is.matrix(x) ) ) )
  if (any(i)) stop('inputs must be vectors or matrices!')
  
  A = args[[1]]
  if (is.vector(A)) {
    rownamesA = names(A)
    A = diag(A, nrow = length(A), ncol = length(A))
    if (!is.null(rownamesA)) {
      rownames(A) = rownamesA
      colnames(A) = rownamesA
    }
  }
  rownamesA = rownames(A)
  colnamesA = colnames(A)
  
  if (n_args == 1) return( A )
  
  for (k in (2:n_args)) {
    B = args[[k]]
    
    if (is.vector(B)) {
      rownamesB = names(B)
      B = diag(B, nrow = length(B), ncol = length(B))
      if (!is.null(rownamesB)) {
        rownames(B) = rownamesB
        colnames(B) = rownamesB
      }
    }
    rownamesB = rownames(B)
    colnamesB = colnames(B)
    
    if ( !((is.null(rownamesA)) && (is.null(rownamesB))) ) {
      if (is.null(rownamesA)) rownamesA = character(nrow(A))
      if (is.null(rownamesB)) rownamesB = character(nrow(B))
      rownamesA = c(rownamesA, rownamesB)
    }
    if ( !((is.null(colnamesA)) && (is.null(colnamesB))) ) {
      if (is.null(colnamesA)) colnamesA = character(ncol(A))
      if (is.null(colnamesB)) colnamesB = character(ncol(B))
      colnamesA = c(colnamesA, colnamesB)
    }
    A = rbind( cbind(A, matrix(FALSE, nrow = nrow(A), ncol = ncol(B))),
               cbind(matrix(FALSE, nrow = nrow(B), ncol = ncol(A)), B) )
    
  }
  
  rownames(A) = rownamesA
  colnames(A) = colnamesA
  return(A)
}

#' Block matrices
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' This helper function coerces an array to a (block) matrix.
#'
#' @param x    Vector, matrix or array. Vectors are coerced to one dimensional arrays and matrices
#'             are of course treated as 2-dimensional arrays.
#' @param rows,cols integer vectors. These two vectors define a partition of the "dimensions" \code{(1,...,n)}, where
#'             \code{n} is the number of dimensions of \code{x} (i.e. \code{length(dim(x))}).
#'             If either of the two is missing, then the complement is used. At least one of the arguments "rows" and
#'             "cols" has to be given.
#'
#' @return matrix
#' @export
#'
#' @examples
#' x = 1:4
#' bmatrix(x, rows = 1, cols = integer(0)) # returns an (4,1) matrix
#' bmatrix(x, cols = 1, rows = NULL)       # returns an (1,4) matrix
#'
#' x = test_array(dim = c(2,3))
#' bmatrix(x, cols = 2)    # returns x    (is equivalent to bmatrix(x, rows = 1, cols = 2))
#' bmatrix(x, rows = 2)    # returns t(x) (is equivalent to bmatrix(x, rows = 2, cols = 1))
#' bmatrix(x, rows = integer(0))   # returns an (1,6) matrix
#'
#' x = test_array(dim = c(2,3,4))
#' bmatrix(x, rows = 1)   # is equivalent to: bmatrix(x, rows = 1, cols = c(2,3))
#' bmatrix(x, rows = 1, cols = c(3,2))
#' bmatrix(x, cols = 2)   # is equivalent to: bmatrix(x, rows = c(1,3), cols = 2)
#' bmatrix(x, rows = 1:3) # is equivalent to: bmatrix(x, cols = integer(0))
#' bmatrix(x, rows = c(3,1,2))
#'
#' \dontrun{
#' # the examples below throw an error
#' bmatrix(x, rows = 1, cols = 2)
#' bmatrix(x, rows = c(1,2), cols = c(2,3))
#' bmatrix(x, rows = c(1,2,1), cols = 3)
#' }
bmatrix = function(x, rows = NULL, cols = NULL) {
  # check input parameters
  if (is.null(rows) && is.null(cols)) {
    stop('missing parameters "rows" and "cols"!')
  }
  if (is.vector(x)) {
    x = array(x, dim = length(x))
  }
  if (!is.array(x)) {
    stop('"x" must be a vector, matrix or array!')
  }
  dim.x = dim(x)
  n.dims = length(dim.x)
  dims = (1:n.dims)
  
  # check "rows"
  if (!is.null(rows)) {
    rows = as.integer(as.vector(rows))
    # rows must be subset of dims
    if ( (length(rows) > 0) && ( (min(rows) < 1) || (max(rows) > n.dims) ) ) {
      stop('parameter "rows" does not correspond to a selection of dimensions (of "x")!')
    }
  }
  
  # check "cols"
  if (!is.null(cols)) {
    cols = as.integer(as.vector(cols))
    # cols must be subset of dims
    if ( (length(cols) > 0) && ( (min(cols) < 1) || (max(cols) > n.dims) ) ) {
      stop('parameter "cols" does not correspond to a selection of dimensions (of "x")!')
    }
  }
  
  if (is.null(rows)) {
    if (length(cols) > 0) {
      rows = dims[-cols]
    } else {
      rows = dims
    }
  }
  if (is.null(cols)) {
    if (length(rows) > 0) {
      cols = dims[-rows]
    } else {
      cols = dims
    }
  }
  
  #  cat('rows:',rows,'cols:',cols,'\n')
  
  # check that the union of "rows" and "cols" = "dims"
  if (!isTRUE(all.equal(dims, sort(c(rows,cols))))) {
    stop('the parameters "rows" and "cols" do not correspond to a partition of the set of dimensions (of "x")!')
  }
  
  # coerce x to a matrix
  x = aperm(x, c(rows,cols))
  dim(x) = c(prod(dim.x[rows]),prod(dim.x[cols]))
  
  return(x)
}


#' Block Toeplitz matrix
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Construct a block Toeplitz matrix from two 3-dimensional arrays \code{R} and \code{C}. The array
#' \code{R} determines the first block row and \code{C} the first block column of the result.
#' The \eqn{(i,j)}-th block of the Toeplitz matrix is \code{R[,,j-i+1]} for \eqn{j\geq i}{j \ge i}
#' and \code{C[,,i-j+1]} for \eqn{j < i}. In particular, note that the \eqn{(1,1)} block
#' is set to \code{R[,,1]} (while \code{C[,,1]} is ignored).
#'
#' If only the argument \code{R} is provided then \code{C} is set to \code{C = aperm(R,c(2,1,3)))}. So
#' \code{btoeplitz(R = X)} is equivalent to \code{btoeplitz(R = X, C = aperm(X,c(2,1,3)))}. Analogously
#' \code{btoeplitz(C = X)} is equivalent to \code{btoeplitz(R = aperm(X,c(2,1,3)), C = X)}.
#' In these cases \eqn{p=q} must hold. The so constructed Toeplitz matrix is
#' symmetric if and only if \code{X[,,1]} is symmetric. Note also that even in the case where only \code{C}
#' is supplied the "\code{R} wins" rule holds, i.e. \code{btoeplitz(C = X)} sets the \eqn{(1,1)} block
#' equal to \code{t(X[,,1])}.
#'
#' @param R Array with dimensions \eqn{(p,q,n)}.
#'   Corresponds to the first (block-) row, containing \eqn{n} matrices of dimension \eqn{(p \times q)}{(p x q)}.
#'   A vector is coerced to an \eqn{(1,1,n)} dimensional array and a \eqn{(p,q)} matrix is interpreted as
#'   an \eqn{(p,q,1)} dimensional array.
#' @param C Array with dimensions \eqn{(p,q,m)}.
#'   Corresponds to the first (block-) column, containing \eqn{m} matrices of dimension \eqn{(p \times q)}{(p x q)}.
#'   A vector is coerced to an \eqn{(1,1,m)} dimensional array and a \eqn{(p,q)} matrix is interpreted as
#'   an \eqn{(p,q,1)} dimensional array.
#'
#' @return Matrix of size \eqn{( pm \times qn )}{( pm x qn )}.
#' @export
#'
#' @examples
#' btoeplitz(0:3)
#' btoeplitz(0:3,-(0:3))
#'
#' btoeplitz(R = test_array(dim=c(2,3,3)), C = -test_array(dim = c(2,3,2)))
#' btoeplitz(R = test_array(dim=c(2,2,1)), C = test_array(dim=c(2,2,0)))
#' btoeplitz(R = test_array(dim=c(2,2,3)))
#' btoeplitz(C = -test_array(dim=c(2,2,3)))
#' # create a symmetric matrix
#' X = test_array(dim= c(2,2,3))
#' X[,,1] = (X[,,1] + t(X[,,1]))/2
#' btoeplitz(R = X)
#' btoeplitz(C = X)
#'
#' \dontrun{
#' # the following examples throw an error
#' btoeplitz(R = test_array(dim=c(2,1,3)), C = -test_array(dim = c(2,2,2)))
#' btoeplitz(R = test_array(dim=c(2,1,3)))
#' btoeplitz(C = test_array(dim=c(2,1,3)))
#' }
btoeplitz = function(R, C) {
  
  # Check for correct inputs: Vector, matrix, or array ####
  # Both missing?
  if ((missing(R)) && (missing(C))) {
    stop('At least one parameter "R" or "C" has to be supplied.')
  }
  
  # R correct?
  if (!missing(R)) {
    if (is.vector(R)) {
      R = array(R, dim = c(1, 1, length(R)))
    }
    if (is.matrix(R)) {
      dim(R) = c(nrow(R), ncol(R), 1)
    }
    if ((!is.array(R)) || (length(dim(R)) != 3)) {
      stop('R must be a vector, a matrix or a 3-dim array!')
    }
    dimR = dim(R)
  }
  
  # C correct?
  if (!missing(C)) {
    if (is.vector(C)) {
      C = array(C, dim = c(1, 1, length(C)))
    }
    if (is.matrix(C)) {
      dim(C) = c(nrow(C), ncol(C), 1)
    }
    if ((!is.array(C)) || (length(dim(C)) != 3)) {
      stop('C must be a vector, a matrix or a 3-dim array!')
    }
    dimC = dim(C)
  }
  
  # If one argument is missing, replace it such that we obtain a symmetric Toeplitz matrix ####
  
  if (missing(R)) {
    if (dim(C)[1] != dim(C)[2]) stop('if "R" is missing then dim(C)[1]=dim(C)[2] must hold!')
    # create a 'symmetric' block toeplitz matrix
    R = aperm(C, c(2, 1, 3))
    dimR = dim(R)
  }
  
  if (missing(C)) {
    if (dim(R)[1] != dim(R)[2]) stop('if "C" is missing then dim(R)[1]=dim(R)[2] must hold!')
    # create a 'symmetric' block toeplitz matrix
    C = aperm(R, c(2, 1, 3))
    dimC = dim(C)
  }
  
  # Check for compatible dimensions ####
  if (any(dimR[1:2] != dimC[1:2])) {
    stop('Dimensions of "R" and "C" are not compatible.')
  }
  
  # Define integer valued parameters (dimensions) ####
  n_rows = dimR[1] # p
  n_cols = dimR[2] # q
  n_depth_R = dimR[3] # n
  n_depth_C = dimC[3] # m
  
  # Deal with trivial cases (one array consists of 0 or 1 matrix, i.e. there is no third dimension) ####
  if (n_depth_C == 0)
    return(matrix(0, nrow = 0, ncol = n_cols * n_depth_R))
  if (n_depth_C == 1)
    return(matrix(R, nrow = n_rows, ncol = n_cols * n_depth_R))
  
  if (n_depth_R == 0)
    return(matrix(0, nrow = n_rows * n_depth_C, ncol = 0))
  if (n_depth_R == 1) {
    C[, , 1] = R[, , 1]
    return(matrix(aperm(C, c(1, 3, 2)), nrow = n_rows * n_depth_C, n_cols))
  }
  
  # Concatenate C and R, and put them in the correct order (to fill block Toeplitz matrix in column-major order) ####
  CR = array(0, dim = c(n_rows, n_cols, n_depth_C + n_depth_R - 1))
  CR[, , 1:(n_depth_C - 1)] = C[, , n_depth_C:2] # order needs to be inverted because we fill from top to bottom
  CR[, , n_depth_C:(n_depth_C + n_depth_R - 1)] = R
  
  # Create indices
  mat_tmp = matrix(raw(), nrow = n_depth_C, ncol = n_depth_R)
  idx_mat = col(mat_tmp) - row(mat_tmp) + n_depth_C
  
  # Create Block Toeplitz matrix
  T = CR[, , idx_mat, drop = FALSE]
  dim(T) = c(n_rows, n_cols, n_depth_C, n_depth_R)
  T = aperm(T, c(1, 3, 2, 4))
  dim(T) = c(n_depth_C * n_rows, n_depth_R * n_cols)
  return(T)
}


#' Block Hankel matrix
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Construct a block Hankel matrix with (block) entries from a 3-dimensional array R.
#' The \eqn{(i,j)}-th block of the Hankel matrix is \code{R[,,i+j-1]}.
#'
#' @param R 3-dimensional array, matrix or vector. A vector of length \eqn{k} is coerced to a
#'          \eqn{(1,1,k)}-dimensional array and a \eqn{(p,q)} matrix is treated as an array of
#'          dimension \eqn{(p,q,1)}.
#' @param d determines the number of block rows and columns.
#'   Suppose that \code{R} is an array of size \eqn{(p,q,k)}.
#'   If \eqn{d=(m,n)} then \code{bhankel} returns a block Hankel matrix with \eqn{m} block rows and \eqn{n} block columns.
#'   If \eqn{d=m} then a Hankel matrix with \eqn{m} block rows and
#'   \eqn{n=\max(k+1-m,1)}{n=max(k+1-m,1)} block columns is returned.
#'   In the default case \code{d = NULL} the number of block rows is \eqn{m=(k+1)/2} for odd \eqn{k}
#'   and \eqn{m=(k/2+1)} else and the number of block columns is set to \eqn{n=\max(k+1-m,1)}{n=max(k+1-m,1)}.
#'   In the case \eqn{(m+n-1)>k} the array \code{R} is padded with zeros.
#'
#' @return Matrix of size \eqn{(pm\times qn)}{(pm x qn)}.
#' @export
#'
#' @examples
#' bhankel(1:6)
#' bhankel(1:6, d = 3)
#' bhankel(letters[1:6], d = c(3,7))  # note the "zero padding"
#' bhankel(test_array(dim = c(2,2,6)))
#' bhankel(test_array(dim = c(3,2,6)), d = 3)
#' bhankel(test_array(dim = c(2,2,6)), d = c(3,7))
#' bhankel(test_array(dim = c(1,2,6)), d = c(1,2))
#' bhankel(test_array(dim = c(2,3)), d = c(2,2))
#' bhankel(test_array(dim = c(2,2,6)), d = c(3,0))
#' bhankel(test_array(dim = c(2,2,0)), d = c(2,2))
#'
#' \dontrun{
#' # the following examples throw an error
#' bhankel(1:5, d = c(-1,1))
#' bhankel(test_array(dim = c(2,3,2,1)))
#' }
bhankel = function(R, d = NULL) {
  if (is.vector(R)) R = array(R,dim=c(1,1,length(R)))
  if (is.matrix(R)) R = array(R,dim=c(nrow(R),ncol(R),1))
  dim = dim(R)
  if (length(dim)!=3) stop('"R" must be a vector, matrix or 3-dimensional array!')
  p = dim[1]
  q = dim[2]
  k = dim[3]
  # if (k==0) stop('"R" is an empty array!')
  
  if ( is.null(d) ) {
    d = ceiling((k+1)/2)
  }
  d = as.integer(as.vector(d))
  if (min(d)<0) stop('"d" has negative entries!')
  if (length(d)==1) d = c(d,max(k+1-d,1))
  d = d[1:2]
  m = d[1]
  n = d[2]
  
  if (min(c(m,n,p,q))==0) return(matrix(0,nrow = m*p,ncol=n*q))
  
  # pad with zeros
  if (k < (m+n-1)) R = array(c(R,double((m+n-1-k)*p*q)),dim=c(p,q,m+n-1))
  
  if (m==1) return(matrix(R[,,1:n],nrow=p,ncol=q*n))
  
  if (n==1) {
    return(matrix(aperm(R[,,1:m,drop=FALSE],c(1,3,2)),nrow=m*p,q))
  }
  
  # j+i-1
  ji = matrix(1:n,nrow=m,ncol=n,byrow = TRUE) + matrix(0:(m-1),nrow=m,ncol=n)
  T = array(0,dim=c(p,q,m*n))
  T = R[,,ji,drop=FALSE]
  dim(T) = c(p,q,m,n)
  T = aperm(T,c(1,3,2,4))
  dim(T) = c(m*p,n*q)
  return(T)
}



# Linear Indices vs Matrix indices ##############################################

#' Transform between Linear Index and Matrix Indices
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' These functions are as in MATLAB.
#' \code{ind2sub()} transforms a linear index to a row and column index for a matrix of given size.
#' \code{sub2ind()} transforms a matrix index, (row, col) to a a linear index (in terms of columns).
#'
#' @param dim Integer vector of size 2. Matrix dimensions.
#' @param ind Integer. Linear index.
#' @param row Integer. Row index.
#' @param col Integer. Column index.
#'
#' @return \code{ind2sub()} returns for given linear index, a matrix index (row, col).
#'   \code{sub2ind()} returns for a given matrix index (row, col), a linear index (column-major).
#'
#' @export
#'
#' @examples
#' A = matrix(1:(3*4), 3, 4)
#' A
#'
#' ind2sub(c(3,4), 7)
#' sub2ind(c(3,4), 2, 3)
#' @name idx_trafo
ind2sub = function(dim, ind){
  row = ((ind-1) %% dim[1]) + 1
  col = floor((ind-1) / dim[1]) + 1
  return(c(row, col))
}


#' @export
#' @rdname idx_trafo
sub2ind = function(dim, row, col){
  ind = (col-1)*dim[1] + row
  return(ind)
}

# QL and LQ decomposition #### 
#' QL and LQ Decomposition
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Returns the QL and LQ decomposition of a matrix with non-negative "diagonal" elements in the QL and LQ decompositions.
#' Only works if no column pivoting occurs in the QR decomposition \code{\link{qr}}.
#'
#' @section Implementation of the QL decomposition using the QR decomposition:
#' MORE DOCU PLEASE
#'
#' @note
#' Shouldn't export this function when publishing package!
#'
#' @param x Matrix.
#' @param ... Other arguments for \code{\link{qr}}
#'
#' @return List with two elements:
#' \item{q}{Orthogonal matrix}
#' \item{l}{Lower triangular matrix with non-negative diagonal elements.}
#' 
#' @rdname ql_decomposition
#' @name QL and LQ decomposition
#' @export
#'
#' @examples
#' # QL Decomposition #############
#' 
#' set.seed(1803)
#' 
#' # Tall matrix
#' x = matrix(stats::rnorm(5*3), 5, 3)
#' out = ql_decomposition(x)
#' all.equal(x, out$q %*% out$l)
#' all.equal(diag(3), t(out$q) %*% out$q)
#' out$l
#'
#' # Wide matrix
#' x = matrix(stats::rnorm(5*3), 3, 5)
#' out = ql_decomposition(x)
#' all.equal(x, out$q %*% out$l)
#' all.equal(diag(3), t(out$q) %*% out$q)
#' out$l
#'
#' # Square matrix
#' x = matrix(stats::rnorm(5*3), 3, 3)
#' out = ql_decomposition(x)
#' all.equal(x, out$q %*% out$l)
#' all.equal(diag(3), t(out$q) %*% out$q)
#' out$l
#'
#' # LQ Decomposition #############
#'
#' # Tall matrix
#' x = matrix(stats::rnorm(5*3), 5, 3)
#' out = lq_decomposition(x)
#' all.equal(x, out$l %*% out$q)
#' all.equal(diag(3), out$q %*% t(out$q))
#' out$l
#'
#' # Wide matrix
#' x = matrix(stats::rnorm(5*3), 3, 5)
#' out = lq_decomposition(x)
#' all.equal(x, out$l %*% out$q)
#' all.equal(diag(3), out$q %*% t(out$q))
#' out$l
#'
#' # Square matrix
#' x = matrix(stats::rnorm(5*3), 3, 3)
#' out = lq_decomposition(x)
#' all.equal(x, out$l %*% out$q)
#' all.equal(diag(3), out$q %*% t(out$q))
#' out$l
#' 
#' # reset seed
#' set.seed(NULL)
ql_decomposition = function(x,...) {
  # dimensions
  m = nrow(x)
  n = ncol(x)
  
  wide = n > m
  
  # Check inputs
  if (min(m,n)<=0) {
    stop('x is an empty matrix!')
  }
  
  x = x[, n:1, drop = FALSE]
  
  qr_x = qr(x,...)
  if (any(qr_x$pivot != 1:n)){
    stop("QL implementation via QR does not work if QR decomposition pivots columns.")
  }
  
  r = qr.R(qr_x)
  q = qr.Q(qr_x)
  
  if (wide) {
    l = r[m:1, n:1, drop = FALSE]
    s = diag(sign(diag(l[1:m, (n-m+1):n, drop = FALSE])))
    l = s %*% l
    q = q[, m:1, drop = FALSE] %*% s
  } else {
    l = r[n:1, n:1, drop = FALSE]
    s = diag(sign(diag(l)))
    l = s%*% l
    q = q[, n:1, drop = FALSE] %*% s
  }
  
  return(list(l = l, q = q))
}



#' @rdname ql_decomposition
#' @export
lq_decomposition = function(x, ...){
  
  # Dimensions
  m = nrow(x)
  n = ncol(x)
  
  wide = n > m
  
  # QR of transpose of x (in order to obtain the LQ of x)
  tx = t(x)
  qr_tx = qr(tx, ...)
  
  if (any(qr_tx$pivot != 1:m)){
    stop("This function does not work if pivoting is used by the base::qr() function.
         Only way to resolve this error is using a version without pivoting
         (of which BF is not aware to exist in R; but it does so in scipy library)")
  }
  
  r = qr.R(qr_tx)
  q = qr.Q(qr_tx)
  
  # wide = t(x) is tall
  if (wide){
    s = diag(sign(diag(r)))
    r = s %*% r
    q = q %*% s
  } else {
    # t(x) is wide and n is the "small" index
    s = diag(sign(diag(r[1:n, 1:n, drop = FALSE])))
    r = s %*% r
    q = q %*% s
  }
  return(
    list(l = t(r),
         q = t(q))
  )
}

# test_ objects (lmfd, rmfd, polm, stsp) #### 

#' Create Test Rational Matrix in LMFD Form
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' This simple tool may be used to create a random, \eqn{(m,n)}-dimensional, rational matrix in 
#' LMFD form 
#' \deqn{x(z)=a^{-1}(z) b(z)}{x(z) = a^{-1}(z) b(z)}
#' The degrees of the polynomials \eqn{a(z), b(z)} is denoted with \eqn{p} and \eqn{q} respectively.
#' 
#' We require \eqn{m>0} and \eqn{p\geq 0}{p\ge 0}. The left factor \eqn{a(z)} 
#' is normalized as \eqn{a(0)=I_m}{a(0)=I} where \eqn{I_m}{I} denotes 
#' the \eqn{(m,m)}-dimensional identity matrix.  
#' 
#' The user may prescribe lower bounds for the moduli of the zeroes and/or poles of the rational matrix.
#' In this case the procedure simply generates (up to n.trials) random matrices until a matrix is found 
#' which satisfies the constraint. The standard deviation of the normal distribution, which is used to 
#' generate the random entries, is decreased in each step. Of course this is a very crude method and 
#' it may fail or need a very large number of randomly generated matrices.
#'
#' @param dim integer vector \code{c(m,n)}. 
#' @param degrees integer vector \code{c(p,q)}.
#' @param digits integer, if non NULL then the randomly generated numbers are rounded to 
#'               "digits" number of decimal places.
#' @param bpoles lower bound for the moduli of the poles of the rational matrix (or NULL).
#' @param bzeroes lower bound for the moduli of the zeroes of the rational matrix (or NULL). 
#'                This parameter is ignored for non-square matrices (m != n).
#' @param n.trials maximum number of trials.
#'
#' @return \code{\link{lmfd}} object, which represents the generated rational matrix \eqn{x(z)} in 
#'         LMFD form. 
#' @export
#' 
#' @examples 
#' ### generate a random (2 x 2) rational matrix in LMFD form with degrees p=1 and q =1
#' ### we require that the matrix has no poles and no zeroes within the unit circle!
#' x = try(test_lmfd(dim = c(2,2), degrees = c(1,1), digits = 2, bpoles = 1, bzeroes = 1))
#' if (!inherits(x, 'try-error')) {
#'    print(x)
#'    print(abs(poles(x)))
#'    print(abs(zeroes(x)))
#' }
test_lmfd = function(dim = c(1,1), degrees = c(1,1), digits = NULL, 
                     bpoles = NULL, bzeroes = NULL, n.trials = 100) {
  # check input parameter "dim"
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (dim[1] <= 0) || (dim[2] < 0)) {
    stop('argument "dim" must be a vector of integers with length 2, dim[1] > 0 and dim[2] >= 0!')
  }
  # check input parameter "degree"
  degrees = as.integer(degrees) # note: as.integer converts to vector!
  if ((length(degrees) != 2) || (degrees[1] < 0) || (degrees[2] < (-1))) {
    stop('argument "degrees" must be a vector of integers with length 2, degrees[1] >= 0 and degrees[2] >= -1!')
  }
  
  m = dim[1]
  n = dim[2]
  p = degrees[1]
  q = degrees[2]
  
  if (p == 0) bpoles = NULL
  if ( (m != n) || (q == 0) ) bzeroes = NULL
  
  i.trial = 0 
  err = TRUE
  sd = 1
  while ( (i.trial < n.trials) && (err) ) {
    a = cbind(diag(m), matrix(stats::rnorm(m*m*p, sd = sd), nrow = m, ncol = m*p))
    dim(a) = c(m,m,p+1)
    b = matrix(stats::rnorm(m*n*(q+1), sd = sd), nrow = m, ncol = n*(q+1))
    dim(b) = c(m,n,q+1)
    if (!is.null(digits)) {
      a = round(a, digits)
      b = round(b, digits)
    }
    x = lmfd(a,b)
    
    err = FALSE
    if ( !is.null(bpoles) ) {
      err = try(min(abs(poles(x, print_message = FALSE))) <= bpoles, silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    if ( (!err) && (!is.null(bzeroes)) ) {
      err = try((min(abs(zeroes(x, print_message = FALSE))) <= bzeroes), silent = TRUE)
      if (inherits(err, 'try-error')) err = TRUE
    }
    i.trial = i.trial + 1
    sd = sd/1.1
  }
  if (err) {
    stop('Could not generate a suitable rational matrix with ', n.trials, ' trials!')
  }
  return(x)
}

# obsolote / old version 
# #' @rdname test_array
# #' @export
# test_polm = function(dim, random = FALSE){
#   # check input parameter "dim"
#   dim = as.vector(as.integer(dim))
#   if ((length(dim) != 3)) {
#     stop('"dim" must be a vector of integers with length 3!')
#   }
#   dim[3] = dim[3] + 1
#   if ( min(dim) < 0) {
#     stop('"dim[1], dim[2] and dim[3]+1" must be non negative!')
#   }
#   
#   if (dim[1]*dim[2]*dim[3] == 0) {
#     return(polm(array(0, dim = dim)))
#   }
#   if (random){
#     polm(array(stats::rnorm(prod(dim)), dim = dim))
#   } else {
#     polm(test_array(dim) - 1)
#   }
# }


# only real matrices
# col_end_matrix changes the degree
# col_end_matrix for degree -1 columns is ignored

#' Create Test Polynomial Matrix
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' This simple tool creates (random) polynomial matrices for testing purposes. 
#' 
#' Note that the desired parameters \code{value_at_0} and \code{column_end_matrix} 
#' may be in conflict to the desired degree(s). See the examples below.
#' 
#' The matrices \code{column_end_matrix} and \code{value_at_0} may contain \code{NA}'s. 
#' These entries are replaced by randomly generated numbers. 
#' 
#' The user may prescribe lower bounds for the moduli of the zeroes 
#' of the polynomial. In this case the procedure simply generates (up to \code{n.trials}) 
#' random matrices until a matrix is found which satisfies the constraint. The standard deviation 
#' of the normal distribution, which is used to generate the random entries, is decreased in each step. 
#' Of course this is a very crude method and it may fail or need a very large number of randomly 
#' generated matrices.
#'
#' @inheritParams test_lmfd
#' @param dim two dimensional vector of non negative integers, determines the dimension of 
#'            the polynomial matrix to be created. If the prescribed number of rows or number of columns 
#'            is zero then an "empty" polynomial matrix is generated. In this case all 
#'            parameters below are ignored. 
#' @param degree desired degree(s), either a scalar, a vector or a matrix. In the matrix case 
#'               \code{degree} prescribes the degrees of the entries of the polynomial matrix. 
#'               A vector \code{degree} defines the column degrees of the matrix, i.e. 
#'               the respective maximal degrees within the columns, and a scalar \code{degree} 
#'               determines the maximum degree of all entries of the matrix. Of course 
#'               \code{degree} has to be compatible with the parameter \code{dim}. 
#'               If the desired degree is \eqn{-1} then a zero polynomial is generated and 
#'               the parameters below are ignored. 
#' @param random If TRUE the coefficents are generated by drawing from a normal distribution.
#'               If FALSE then the coefficient of the \eqn{k}-th power \eqn{z^k} of 
#'               the (i,j)-th entry is set equal to "\eqn{ijk}". In this case the 
#'               parameters below are ignored!
#' @param col_end_matrix desired column end matrix (or \code{NULL}).  
#' @param value_at_0 desired value of the polynomial at \eqn{z=0} (or \code{NULL}).
#'
#' @return \code{\link{polm}} object, which represents the generated polynomial matrix. 
#' @export
#'
#' @examples
#' ### "empty" polynomials, number of rows or number of columns is zero. 
#' test_polm(dim = c(0,0))
#' test_polm(dim = c(0,2))
#' test_polm(dim = c(3,0))
#' 
#' ### (3,3) polynomial of degree -1 (i.e. a(z)=0) 
#' test_polm(dim = c(3,3), degree = -1)
#' 
#' ### (3,3) polynomial with degree 1 
#' test_polm(dim = c(3,3), degree = 1) %>% print(format = 'c')
#' 
#' ### (3,3) polynomial with column degrees -1, 0 and 1 
#' test_polm(dim = c(3,3), degree = c(-1,0,1))  %>% print(format = 'c')
#' 
#' ### random, (3,3) polynomial with prescribed degrees for each element 
#' deg = matrix(c(0, 1, 2,
#'                1,-1,-1,
#'                1, 1, 1), nrow = 3, ncol = 3, byrow = TRUE)
#' print(deg)
#' 
#' a = test_polm(dim = c(3,3), degree = deg, random = TRUE)
#' print(a, digits = 2, format = 'c')
#' 
#' ### random, (3,3) polynomial with prescribed column degree and column end matrix
#' cm = matrix(NA_real_, nrow = 3, ncol = 3)
#' cm[lower.tri(cm, diag = FALSE)] = 0
#' a = test_polm(dim = c(3,3), degree = c(0,1,2), random = TRUE, 
#'               digits = 2, col_end_matrix = cm)
#' print(a, digits = 2, format = 'c')
#' print(degree(a, which = 'column'))
#' print(col_end_matrix(a))
#' 
#' ### the parameters column_end_matrix and value_at_zero 
#' ### may be in conflict with the prescribed degree(s). 
#' 
#' # E.g. if we set the second column of "cm" equal to zero
#' cm[, 2] = 0
#' a = test_polm(dim = c(3,3), degree = c(0,1,2), random = TRUE, 
#'               digits = 2, col_end_matrix = cm)
#' print(a, digits = 2, format = 'c')
#' 
#' # then the generated polynomial has column degrees 0,0,2 
#' # and the column end matrix is not upper triangular!
#' print(degree(a, which = 'column'))
#' print(col_end_matrix(a))
#' 
#' ### here we set a(0) equal to the identity matrix and 
#' ### require that a(z) has no zeroes within the unit circle
#' a = try(test_polm(dim = c(3,3), degree = 2, random = TRUE, 
#'         digits = 2, value_at_0 = diag(3), bzeroes = 1))
#' if (!inherits(a, 'try-error')) {
#'   print(a, digits = 2, format = 'c')
#'   print(abs(zeroes(a)))
#' }

test_polm = function(dim = c(1,1), degree = 0, random = FALSE, digits = NULL, col_end_matrix = NULL, 
                     value_at_0 = NULL, bzeroes = NULL, n.trials = 100) {
  # check input parameter "dim"
  dim = as.integer(dim) # note: as.integer converts to vector!
  if ((length(dim) != 2) || (min(dim) < 0)) {
    stop('argument "dim" must be a vector of non negative integers with length 2!')
  }
  
  # if an empty matrix is desired, ignore all other parameters
  if (prod(dim) == 0) {
    return(polm(array(0, dim = c(dim,0))))
  }
  
  # check input parameter "degree"
  if (is.vector(degree)) {
    degree = as.integer(degree) 
    if (length(degree) == 1) degree = rep(degree, dim[2])
    if (length(degree) != dim[2]) stop('parameters "degree" and "dim" are not compatible.')
    degree = matrix(degree, nrow = dim[1], ncol = dim[2], byrow = TRUE)
  }
  if (is.matrix(degree)) {
    # coerce to integer matrix. note that as.integer() returns a vector!
    degree = matrix(as.integer(degree), nrow = nrow(degree), ncol = ncol(degree))
  }
  if ( (!is.integer(degree)) || (!is.matrix(degree)) || 
       any(dim(degree) != dim) || any(is.na(degree)) || any(degree < -1) ) {
    stop('argument "degree" must be a scalar, vector or matrix of integers (>= -1), compatible to "dim".')
  }
  
  p = max(degree)
  if (p == (-1)) {
    # zero polynomial
    return(polm(array(0, dim = c(dim,0))))
  }
  # compute column degrees
  col_degree = apply(degree, MARGIN = 2, FUN = max)
  
  # create a polynomial with 'fixed' coefficients 
  if (!random) {
    x = test_array(dim = c(dim, p+1)) - 1
    # impose the desired degrees!
    for (i in (1:dim[1])) {
      for (j in (1:dim[2])) {
        if (degree[i,j] < p) {
          x[i,j,(degree[i,j]+2):(p+1)] = 0
        }
      }
    }
    return(polm(x))
  }
  
  # create a random polynomial   
  
  # create an array with NA's for the "free" parameters
  x0 = array(NA_real_, dim = c(dim, p+1))
  
  # impose the desired degrees!
  for (i in (1:dim[1])) {
    for (j in (1:dim[2])) {
      if (degree[i,j] < p) {
        x0[i, j, (degree[i,j]+2):(p+1)] = 0
      }
    }
  }
  
  # impose column end matrix 
  if (!is.null(col_end_matrix)) {
    # check input parameter col_end_matrix 
    if ( !is.numeric(col_end_matrix) || !is.matrix(col_end_matrix) || 
         any(dim(col_end_matrix) != dim) ) {
      stop('argument "col_end_matrix"  is not compatible')
    }
    if (min(col_degree) < 0) { 
      stop('some column degrees are negative, but column end matrix has been given')
    }
    for (i in (1:dim[2])) {
      x0[,i,col_degree[i]+1] = col_end_matrix[,i]
    }
  }
  
  # impose value at z=0
  if (!is.null(value_at_0)) {
    # check input parameter col_end_matrix 
    if ( !is.numeric(value_at_0) || !is.matrix(value_at_0) || 
         any(dim(value_at_0) != dim) ) {
      stop('argument "value_at_0"  is not compatible')
    }
    if (min(col_degree) < 0) { 
      stop('some column degrees are negative, but value at z=0 has been given')
    }
    x0[,,1] = value_at_0
  }
  
  i = is.na(x0)
  n_theta = sum(i) # number of "free" parameters 
  
  # for non-square polynomials, or polynoials with degree p=0, ignore the parameter "bzeroes"
  if ( (dim[1] != dim[2]) || (p <= 0) ) bzeroes = NULL
  
  if (n_theta == 0) {
    x = polm(x0)
    if (!is.null(bzeroes)) {
      err = try(min(abs(zeroes(x, print_message = FALSE))) <= bzeroes, silent = TRUE)
      if ( inherits(err, 'try-error') ) err = TRUE
      if (err) {
        stop('the zeroes of the generated polynomial do not satisfy the constraint (bzeroes)!')
      }
    }
    return(x)
  }
  
  i.trial = 0 
  err = TRUE
  sd = 1
  while ( (i.trial < n.trials) && (err) ) {
    theta = stats::rnorm(n_theta, sd = sd)
    if (!is.null(digits)) theta = round(theta, digits)
    x0[i] = theta
    x = polm(x0)
    err = FALSE
    
    if ( !is.null(bzeroes) ) {
      err = try(min(abs(zeroes(x, print_message = FALSE))) <= bzeroes, silent = TRUE)
      if ( inherits(err, 'try-error') ) err = TRUE
    }
    sd = sd/1.1
    i.trial = i.trial + 1
  }
  
  if (err) {
    stop('Could not generate a suitable rational matrix with ', n.trials, ' trials!')
  }
  return(x)
}


#' Match Two Vectors
#'
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Given two vectors \code{x,y} of the same length, the routine \code{match_vectors} returns two 
#' permutations \code{i,j} such that \code{x[i]} matches \code{y[j]} as best as possible. The 
#' procedure uses a simple, greedy search strategy. In the first step \code{i[1], j[1]} are chosen such 
#' that \code{abs(x[i[1]]-y[j[1]])} is equal to the minimum distance between the entries of \code{x} and \code{y}. 
#' Then the procedure iterates this search strategy with the remaining entries of \code{x} and \code{y}. 
#' 
#' @param x,y two vectors of the same length, \eqn{p} say. 
#'
#' @return A list with components 
#' \item{i,j}{two permutations, i.e. integer vectors of length \code{p} with unique entries.}
#' \item{dist}{numeric vector, \code{dist = abs(x[i]-y[j])}. By construction the entries 
#'             of \code{d} are sorted.}
#' \item{d}{numeric matrix with entries \code{d[k,l] = abs(x[k] - y[l])}.}
#' 
#' @export
#' @keywords internal
#'
#' @examples
#' # Match the roots of two polynomials a1 and a2
#' p = 5
#' a1 = rnorm(p+1)
#' a2 = a1 + rnorm(p+1)*(1e-6) # a2 is a "noisy" copy of a1
#' z1 = polyroot(a1)
#' z2 = polyroot(a2)
#' out = match_vectors(z1, z2)
#' print(cbind(z1[out$i], z2[out$j], out$dist))
#' 
#' # A polynomial with real coefficients has pairs of complex conjugate roots
#' # However, the roots returned by "polyroot" in general do not have this property!
#' # Match the roots and their complex conjugates
#' out = match_vectors(z1, Conj(z1))
#' print(cbind(out$i, out$j, out$dist, 
#'             Re(z1[out$i]), Re(Conj(z1)[out$j]), 
#'             Im(z1[out$i]), Im(Conj(z1)[out$j])))
#' 
match_vectors = function(x, y = Conj(x)) {
  x = as.vector(x)
  y = as.vector(y)
  p = length(x)
  if (length(y) != p) stop('non compatible vectors "x" and "y"!')
  
  match = matrix(integer(2*p), nrow = p, ncol = 2)
  dist = numeric(p)
  
  if (p == 0) return(list(i = integer(0), j = integer(0), dist = dist, d = matrix(0, nrow = 0, ncol = 0)))
  
  # matrix of distances abs(x[i] - y[j])
  d = abs(matrix(x, nrow = p, ncol = p) - matrix(y, nrow = p, ncol = p, byrow = TRUE))
  
  if (p == 1) return(list(i = 1L, j = 1L, dist = d[1,1], d = d))
  
  d_k = d
  i = (1:p)
  j = (1:p)
  
  for (k in (1:p)) {
    i0 = which.min(apply(d_k, MARGIN = 1, FUN = min))
    j0 = which.min(d_k[i0,])
    match[k,] = c(i[i0], j[j0])
    dist[k] = d_k[i0,j0]
    # skip i0,j0
    if (k < p) {
      d_k = d_k[-i0,,drop = FALSE]
      d_k = d_k[,-j0,drop = FALSE]
      i = i[-i0]
      j = j[-j0]
    }
  }
  return(list(i = match[,1], j = match[,2], dist = dist, d = d))  
}


