# arithmetics for rational matrices ##############################################  

# internal function 
# upgrade objects to common class 
# * numeric matrices are coereced to polm objects 
# * lmfd objcts are coerced to stsp objects
# * pseries to zvalues throws an error!
# force = TRUE means coercion even if only one argument is supplied
upgrade_objects = function(force = TRUE, ...) {
  args = list(...)
  n_args = length(args)
  if (n_args == 0) {
    # return empty list
    return(as.vector(integer(0), mode = 'list')) 
    }
  
  # skip NULL's
  not_null = sapply(args, FUN = function(x) {!is.null(x)})
  args = args[not_null]
  n_args = length(args)
  if (n_args == 0) {
    # return empty list
    return(as.vector(integer(0), mode = 'list'))
  }
  if ((n_args == 1) && (!force)) {
    # just return arg(s)
    return(args)
  }
  
  classes = sapply(args, FUN = function(x) {class(x)[1]} )
  grades = match(classes, c('polm','lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)
  
  # coerce arguments to a common class = "highest" class
  k = which.max(grades)
  max_grade = grades[k]

  # coerce matrix to polm object  
  if (max_grade == 0) {
    # this should not happen ?!
    obj = try(polm(args[[k]]))
    if (inherits(obj, 'try-error')) stop('could not coerce object to "polm" object')
    args[[k]] = obj
    max_grade = 1
    classes[k] = 'polm'
    grades[k] = 1
  }
  
  # coerce lmfd or rmfd object to stsp object  
  if (max_grade %in% c(2,3)) {
    obj = as.stsp(args[[k]])
    args[[k]] = obj
    max_grade = 4
    classes[k] = 'stsp'
    grades[k] = 4
  }
  
  for (i in (1:n_args)) {
    if ((grades[i] < max_grade) && ( i != k )) {
      if (grades[i] == 0) {
        obj = try(polm(args[[i]]))
        if (inherits(obj, 'try-error')) stop('could not coerce object to "polm" object')
        args[[i]] = obj
      }
      
      if (max_grade == 4) {
        obj = try(as.stsp(args[[i]]))
        if (inherits(obj, 'try-error')) stop('could not coerce object to "stsp" object')
        args[[i]] = obj
      }
      
      if (max_grade == 5) {
        obj = pseries(args[[i]], lag.max = dim(args[[k]])[3])
        args[[i]] = obj
      }
      
      if (max_grade == 6) {
        if (grades[i] == 5) stop('cannot coerce power series to frequency response')
        z = attr(args[[k]], 'z')
        obj = zvalues(args[[i]], z = z)
        args[[i]] = obj
      }
    }
  }
  
  return(args)
}




# internal function 
# create a diagonal rational matrix from a scalar 
# brauche ich das noch ????
diag_object = function(e, d) {
  cl = class(e)[1]

  if (cl == 'polm') {
    e = unclass(prune(e, tol = 0))
    e = as.vector(e)
    if ((length(e) == 0) || (d == 0)) {
      return(polm(matrix(0, nrow = d, ncol = d)))
    }
    a = matrix(0, nrow = length(e), ncol = d*d)
    a[, diag(matrix(1:(d*d), nrow = d, ncol = d))] = e
    a = t(a)
    dim(a) = c(d,d,length(e))
    return(polm(a))
  }
  
  if (cl == 'stsp') {
    if (d == 0) {
      return(stsp(A = matrix(0, nrow = 0, ncol = 0), 
                  B = matrix(0, nrow = 0, ncol = 0),
                  C = matrix(0, nrow = 0, ncol = 0),
                  D = matrix(0, nrow = 0, ncol = 0)))
    }
    A = e$A
    B = e$B
    C = e$C
    D = e$D
    A = do.call(bdiag, args = rep(list(A), d))
    B = do.call(bdiag, args = rep(list(B), d))
    C = do.call(bdiag, args = rep(list(C), d))
    D = diag(x = D[1,1], nrow = d, ncol = d)
    return(stsp(A, B, C, D))
  }
  
  if (cl == 'zvalues') {
    z = attr(e, 'z')
    if ((length(z) == 0) || (d == 0)) {
      a = array(0, dim = c(d,d,length(z)))
      a = structure(a, z = z, class = c('zvalues', 'ratm'))
      return(a)
    }
    e = unclass(e)
    e = as.vector(e)
    a = matrix(0, nrow = length(z), ncol = d*d)
    a[, diag(matrix(1:(d*d), nrow = d, ncol = d))] = e
    a = t(a)
    dim(a) = c(d,d,length(z))
    a = structure(a, z = z, class = c('zvalues', 'ratm'))
    return(a)
  }
  
  stop('unsupported class "', cl, '" for diagonalization of scalars')
}


# internal function
# inflate a scalar to an (m x n) matrix with identical entries 
inflate_object = function(e, m, n) {
  cl = class(e)[1]
  
  if (cl == 'polm') {
    e = unclass(prune(e, tol = 0))
    e = as.vector(e)
    if ((length(e)*m*n) == 0) {
      return(polm(matrix(0, nrow = m, ncol = n)))
    }
    a = array(e, dim = c(length(e), m, n))
    a = aperm(a, c(2,3,1))
    return(polm(a))
  }
  
  if (cl == 'stsp') {
    if ((m*n) == 0) {
      return(stsp(A = matrix(0, nrow = 0, ncol = 0), 
                  B = matrix(0, nrow = 0, ncol = n),
                  C = matrix(0, nrow = m, ncol = 0),
                  D = matrix(0, nrow = m, ncol = n)))
    }
    A = e$A
    s = ncol(A)
    B = e$B
    C = e$C
    D = e$D[1,1]
    if (n < m) {
      A = do.call(bdiag, args = rep(list(A), n))
      B = do.call(bdiag, args = rep(list(B), n))
      C = as.vector(C)
      C = do.call(cbind, args = rep(list(matrix(C, nrow = m, ncol = length(C), byrow = TRUE)), n))
    } else {
      A = do.call(bdiag, args = rep(list(A), m))
      C = do.call(bdiag, args = rep(list(C), m))
      B = as.vector(B)
      B = do.call(rbind, args = rep(list(matrix(B, nrow = length(B), ncol = n)), m))
    }
    D = matrix(D, nrow = m, ncol = n)
    return(stsp(A, B, C, D))
  }
  
  if (cl == 'pseries') {
    max.lag = unname(dim(e)[3])
    if (((max.lag + 1)*m*n) == 0) {
      a = array(0, dim = c(m,n,max.lag+1))
      a = structure(a, class = c('pseries', 'ratm'))
      return(a)
    }
    e = unclass(e)
    e = as.vector(e)
    a = array(e, dim = c(max.lag+1, m, n))
    a = aperm(a, c(2,3,1))
    a = structure(a, class = c('pseries', 'ratm'))
    return(a)
  }
  
  if (cl == 'zvalues') {
    z = attr(e, 'z')
    if ((length(z)*m*n) == 0) {
      a = array(0, dim = c(m,n,length(z)))
      a = structure(a, z = z, class = c('zvalues', 'ratm'))
      return(a)
    }
    e = unclass(e)
    e = as.vector(e)
    a = array(e, dim = c(length(z), m, n))
    a = aperm(a, c(2,3,1))
    a = structure(a, z = z, class = c('zvalues', 'ratm'))
    return(a)
  }

  stop('unsupported class "', cl, '" for inflation of scalars')
}

# internal function
# multiplication of matrix polynomials
# this function performs only basic checks on the inputs!
mmult_poly = function(a, b) {
  # a,b must be two compatible arrays
  da = dim(a)
  db = dim(b)
  if (da[2] != db[1]) stop('arguments are not compatible')
  # if any of the arguments is an empty polynomial, or a polynomial of degree -1
  if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], db[2], 0)))
  
  # skip zero leading coefficients
  if (da[3] > 0) {
    a = a[ , , rev(cumprod(rev(apply(a == 0, MARGIN = 3, FUN = all)))) == 0, drop = FALSE]
    da = dim(a)
  }
  # skip zero leading coefficients
  if (db[3] > 0) {
    b = b[ , , rev(cumprod(rev(apply(b == 0, MARGIN = 3, FUN = all)))) == 0, drop = FALSE]
    db = dim(b)
  }
  # if any of the arguments is an empty polynomial, or a polynomial of degree -1
  if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], db[2], 0)))
  
  pa = da[3] - 1
  pb = db[3] - 1
  
  x = array(0, dim = c(da[1], db[2], pa + pb + 1))
  # the 'convolution' of the coefficients is computed via a double loop
  # of course this could be implemented more efficiently!
  for (i in (0:(pa+pb))) {
    for (k in iseq(max(0, i - pb), min(pa, i))) {
      x[,,i+1] = x[,,i+1] + 
        matrix(a[,,k+1], nrow = da[1], ncol = da[2]) %*% matrix(b[,,i-k+1], nrow = db[1], ncol = db[2])
    }
  }
  return(x)
}


# internal function
# multiplication of two impulse response functions
# this function performs only basic checks on the inputs!
# almost equal to mmult_poly
mmult_pseries = function(a, b) {
  # a,b must be two compatible arrays
  da = dim(a)
  db = dim(b)
  if (da[2] != db[1]) stop('arguments are not compatible')
  # if any of the arguments is an empty pseries, or a pseries with lag.max = -1
  if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], db[2], min(da[3], db[3]))))
  
  lag.max = min(da[3], db[3]) - 1
  # truncate to the minimum lag.max
  # a = a[ , , 1:(lag.max+1), drop = FALSE]
  # b = b[ , , 1:(lag.max+1), drop = FALSE]
  
  x = array(0, dim = c(da[1], db[2], lag.max + 1))
  # the 'convolution' of the impulse response coefficients is computed via a double loop
  # of course this could be implemented more efficiently!
  for (i in (0:lag.max)) {
    for (k in (0:i)) {
      x[,,i+1] = x[,,i+1] + 
        matrix(a[,,k+1], nrow = da[1], ncol = da[2]) %*% matrix(b[,,i-k+1], nrow = db[1], ncol = db[2])
    }
  }
  return(x)
}



# internal function
# elementwise multiplication of matrix polynomials
# this function performs only basic checks on the inputs!
emult_poly = function(a, b) {
  # a,b must be two compatible arrays
  da = dim(a)
  db = dim(b)
  if ( (da[1] != db[1]) || (da[2] != db[2]) ) stop('arguments are not compatible')
  # if any of the arguments is an empty polynomial, or a polynomial of degree -1
  if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], da[2], 0)))
  
  # skip zero leading coefficients
  if (da[3] > 0) {
    a = a[ , , rev(cumprod(rev(apply(a == 0, MARGIN = 3, FUN = all)))) == 0, drop = FALSE]
    da = dim(a)
  }
  # skip zero leading coefficients
  if (db[3] > 0) {
    b = b[ , , rev(cumprod(rev(apply(b == 0, MARGIN = 3, FUN = all)))) == 0, drop = FALSE]
    db = dim(b)
  }
  # if any of the arguments is an empty polynomial, or a polynomial of degree -1
  if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], da[2], 0)))
  
  pa = da[3] - 1
  pb = db[3] - 1
  
  x = array(0, dim = c(da[1], db[2], pa + pb + 1))
  # the 'convolution' of the coefficients is computed via a double loop
  # of course this could be implemented more efficiently!
  for (i in (0:(pa+pb))) {
    for (k in iseq(max(0, i - pb), min(pa, i))) {
      x[,,i+1] = x[,,i+1] + 
        matrix(a[,,k+1], nrow = da[1], ncol = da[2]) * matrix(b[,,i-k+1], nrow = db[1], ncol = db[2])
    }
  }
  return(x)
}


# internal function
# elementwise multiplication of two impulse response functions
# this function performs only basic checks on the inputs!
# almost equal to emult_poly
emult_pseries = function(a, b) {
  # a,b must be two compatible arrays
  da = dim(a)
  db = dim(b)
  if ( (da[1] != db[1]) || (da[2] != db[2]) ) stop('arguments are not compatible')
  # if any of the arguments is an empty pseries, or a pseries with lag.max = -1
  if (min(c(da,db)) == 0) return(array(0, dim = c(da[1], da[2], min(da[3], db[3]))))
  
  lag.max = min(da[3], db[3]) - 1
  # truncate to the minimum lag.max
  # a = a[ , , 1:(lag.max+1), drop = FALSE]
  # b = b[ , , 1:(lag.max+1), drop = FALSE]
  
  x = array(0, dim = c(da[1], db[2], lag.max + 1))
  # the 'convolution' of the impulse response coefficients is computed via a double loop
  # of course this could be implemented more efficiently!
  for (i in (0:lag.max)) {
    for (k in (0:i)) {
      x[,,i+1] = x[,,i+1] + 
        matrix(a[,,k+1], nrow = da[1], ncol = da[2]) * matrix(b[,,i-k+1], nrow = db[1], ncol = db[2])
    }
  }
  return(x)
}


# internal function
# elementwise multiplication of two rational vectors (in stsp form)
emult_stsp_vek = function(a,b) {

  dim_a = dim(a)
  m = dim_a[1]
  s_a = dim_a[3]
  dim_b = dim(b)
  s_b = dim_b[3]
  
  # convert a to a diagonal matrix
  A = a$A
  B = a$B
  C = a$C
  D = a$D
  A = do.call(bdiag, args = rep(list(A), m))
  B = do.call(bdiag, args = rep(list(B), m))
  C = do.call(bdiag, args = lapply(as.vector(1:m, mode = 'list'), 
                                   FUN = function(x) C[x,,drop = FALSE]))
  D = diag(x = as.vector(D), nrow = m, ncol = m)
  da = stsp(A, B, C, D)
  # print(da)
  
  # multiply with b
  ab = da %r% b
  # print(ab)
  # print(pseries(ab) - pseries(a)*pseries(b))

  # controllability matrix  
  Cm = ctr_matrix(ab)
  svd_Cm = svd(Cm, nv = 0)
  # print(svd_Cm$d)
  
  # Cm has rank <= s = s_a+s_b (should we check this?)
  # state transformation
  A = t(svd_Cm$u) %*% ab$A %*% svd_Cm$u
  B = t(svd_Cm$u) %*% ab$B
  C = ab$C %*% svd_Cm$u
  D = ab$D
  s = s_a + s_b

  # print(ctr_matrix(stsp(A,B,C,D)))
  
  # skip the "non-controllable" states 
  ab = stsp(A[1:s,1:s, drop = FALSE], B[1:s,,drop = FALSE], 
            C[,1:s,drop = FALSE], D)
  
  # print(ctr_matrix(ab))
  # print(svd(ctr_matrix(ab))$d)
  
  return(ab)
}


# internal function
# elementwise multiplication of two rational vectors (in stsp form)
emult_stsp = function(a,b) {
  dim_a = dim(a)
  m = dim_a[1]
  n = dim_a[2]
  if ((m*n) == 0) {
    return(stsp(A = matrix(0, nrow = 0, ncol = 0), B = matrix(0, nrow = 0, ncol = n),
                C = matrix(0, nrow = m, ncol = 0), D = matrix(0, nrow = m, ncol = n)))    
  }
  
  s_a = dim_a[3]
  s_b = dim(b)[3]
  
  if ((s_a+s_b)==0) {
    return(stsp(a$A, a$B, a$C, a$D * b$D))
  }
  
  if (n <= m) {
    cols = vector(n, mode = 'list')
    # compute elementwise multiplication of the columns of a and b
    for (i in (1:n)) {
      cols[[i]] = emult_stsp_vek(a[,i], b[,i])
    }
    # print(cols)
    # bind the columns
    ab = do.call(cbind, cols)
    return(ab)
  }
  
  # consider the transposed marices
  ab = t(emult_stsp(t(a), t(b)))
  return(ab)
}

#' Matrix Multiplication of Rational Matrices
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' 
#' @param e1,e2 two rational matrix objects (i.e. \code{\link{polm}}, \code{\link{lmfd}}, \code{\link{rmfd}},  
#'              \code{\link{stsp}}, \code{\link{pseries}} or \code{\link{zvalues}} objects), 
#'              or objects which may be coerced to rational matrix objects. 
#'
#' @return rational matrix object. The class depends on the classes of the arguments \code{e1,e1}.
#' @export
'%r%' = function(e1, e2) {
  
  if ( !( inherits(e1, 'ratm') || inherits(e2, 'ratm') )) {
    stop('one of the arguments must be a rational matrix object (ratm)')
  }
  
  # print(class(e1))
  # print(class(e2))
  
  out = upgrade_objects(force = TRUE, e1, e2)
  e1 = out[[1]]
  e2 = out[[2]]
  
  # print(class(e1))
  # print(class(e2))
  
  d1 = unname(dim(e1))
  cl1 = class(e1)[1]
  gr1 = match(cl1, c('polm','lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)
  d2 = unname(dim(e2))
  cl2 = class(e2)[1]
  gr2 = match(cl1, c('polm','lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)
  
  if (d1[2] != d2[1]) stop('the rational matrices e1, e2 are not compatible (ncol(e1) != nrow(e2))')
  
  # finally do the computations 
  if (cl1 == 'polm') {
    e = polm(mmult_poly(unclass(e1), unclass(e2)))
    e = prune(e, tol = 0)
    return(e)
  }
  
  if (cl1 == 'stsp') {
    A1 = e1$A
    B1 = e1$B
    C1 = e1$C
    D1 = e1$D
    A2 = e2$A
    B2 = e2$B
    C2 = e2$C
    D2 = e2$D
    A = rbind(cbind(A1,                                         B1 %*% C2), 
              cbind(matrix(0,nrow = nrow(A2), ncol = ncol(A1)), A2))
    B = rbind(B1 %*% D2, B2)
    C = cbind(C1, D1 %*% C2)
    D = D1 %*% D2
    # print(A)
    e = stsp(A, B, C, D)
    return(e)
  }
  
  if (cl1 == 'pseries') {
    e = mmult_pseries(unclass(e1), unclass(e2))
    class(e) = c('pseries','ratm')
    return(e)
  }
  
  if (cl1 == 'zvalues') {
    z1 = attr(e1,'z')
    z2 = attr(e2,'z')
    if (!isTRUE(all.equal(z1, z2))) {
      stop('the complex points z of the two frequency responses do not coincide')
    }
    e1 = unclass(e1)
    e2 = unclass(e2)
    e = array(0, dim = c(d1[1], d2[2], length(z1)))
    for (i in (1:length(z1))) {
      e[,,i] = matrix(e1[,,i], nrow = d1[1], ncol = d1[2]) %*% matrix(e2[,,i], nrow = d2[1], ncol = d2[2])
    }
    
    e = structure(e, z = z1, class = c('zvalues', 'ratm'))
    return(e)
  }
  
  stop('this should not happen')
}





#' Arithmetic Ops Group Methods for Rational Matrices
#' 
#' This function was originally part of the R-package \strong{rationalmatrices}.
#' \cr
#' Implements the following basic arithmetic operations on rational matrices
#' \itemize{
#' \item unitary operators \code{'+a'} and \code{'-a'}. 
#' \item the power (\code{'a^k'}) operator for \emph{square, non empty} rational matrices 
#'       \code{a} and integer powers \code{k}.   
#' \item elementwise multiplication (\code{'a * b'}) 
#' \item addition and substraction (\code{'a + b'} and \code{'a - b'})
#' \item elementwise polynomal division (\code{'a \%/\% b'}),
#' \item elementwise polynomial remainder (\code{'a \%\% b'}),
#' }
#' 
#' The unitary operators \code{'+a'} and \code{'-a'} are implemented for all classes. 
#' 
#' The power (\code{'a^k'}) operator is only implemented \emph{square}, rational matrices and integer powers \code{k}.
#' \itemize{
#' \item{\code{a^0} returns the identity matrix (represented by an object of the same class as the input argument \code{a}).} 
#' \item{\code{a^1} simply returns the input arguments \code{a}.} 
#' \item{The case \eqn{k>1} is implemented for all classes. 
#'     However, \code{lmfd} and \code{rmfd} objects are first coerced to \code{stsp} objects.} 
#' \item{The case \eqn{k<0} is implemented for all classes. 
#'     However, the matrix must be non empty and \code{polm} objects and \code{lmfd} and \code{rmfd} objects are first coerced to \code{stsp} objects.}
#' }
#'        
#' For the binary operators `a + b`, `a - b` and `a * b` the two arguments are 
#' first coerced to a common class according to the following scheme 
#' \tabular{rcll}{
#' matrix \tab -> \tab polm \tab (coerce non rational matrices to \code{polm} objects)\cr
#' lmfd   \tab -> \tab stsp \tab (coerce \code{lmfd} objects to \code{stsp} objects)\cr
#' rmfd   \tab -> \tab stsp \tab (coerce \code{rmfd} objects to \code{stsp} objects)\cr
#' polm o stsp  \tab -> \tab stsp \tab \cr
#' polm o pseries  \tab -> \tab pseries \tab \cr
#' polm o zvalues  \tab -> \tab zvalues \tab \cr
#' stsp o pseries  \tab -> \tab pseries \tab \cr
#' stsp o zvalues  \tab -> \tab zvalues \tab \cr
#' }
#' 
#' Note that \code{pseries} objects cannot be (easily) coerced to 
#' \code{zvalues} objects, so these binary operations throw an error, if 
#' one tries to combine an \code{pseries} with a \code{zvalues} object.
#' 
#' If two \code{pseries} objects are combined then they are truncated to the 
#' minimum of the respective number of “lags”. Two \code{zvalues} objects are 
#' only combined if the “z” values are identical. Otherwise an error is thrown.
#' 
#' Of course the two arguments have to be of compatible dimension. 
#' If one of the arguments is a scalar (\eqn{(1,1)} matrix) then this 
#' argument is "expanded" to a compatible matrix with identtical entries. 
#' 
#' Note that the computed statespace realizations are often non minimal! 
#' (This remark also applies for other operations on statespace realizations.)
#' 
#' For the matrix multiplication, see \code{\link{\%r\%}}.
#' 
#' The elementwise polynomal division (\code{'a \%/\% b'}) and the 
#' elementwise polynomial remainder (\code{'a \%\% b'}) are of course only implemented 
#' for polynomial matrices (\code{polm} objects) or (objects which may be coerced to 
#' \code{polm} objects). 
#' 
#' The above remark on scalar arguments also applies for these operations. 
#' 
#' @param e1,e2 At least one of \code{e1, e2} must be an object of 
#' class \code{\link{polm}}, \code{\link{lmfd}}, \code{\link{rmfd}}, \code{\link{stsp}},
#' \code{\link{pseries}} or \code{\link{zvalues}}.
#'
#' @return Rational matrix object.
#' @export
#' 
#' @examples
#' # Multiplication (and division) of a scalar (from left and right)
#' a = test_polm(dim = c(2,2), degree = 3)
#' a
#' -a
#' a * (-1)
#' a %/% 0.5
#'
#' # Addition
#' 2 + a        # 2 is coerced to a constant matrix polynomial
#'
#' # Elementwise remainder
#' a %% 0.5
#' 0.5 %% a
#'
#' # Elementwise division and multiplication with scalar polm
#' z = polm(c(0,1))
#' a %/% z
#' z * a
#' a * z^2
#'
#' # (Non-negative integer) power of univariate polynomial 
#' # (useful for generating a polynomial matrix)
#' z^3
#' matrix(0, nrow = 2) + z * matrix(1, nrow = 2) + z^2 * matrix(2, nrow = 2)
#'
#' # (Non-negative integer) power of quadratic polynomial matrices
#' a^0
#' a^1
#' a^2
#'
Ops.ratm = function(e1, e2) {
  
  # unary operator +/- ###############################################################
  if (missing(e2)) {
    
    if (.Generic == '+') {
      return(e1)
    }
    
    if (.Generic == '-') {
      cl1 = class(e1)[1]
      if (cl1 == 'polm') {
        return(polm(-unclass(e1)))
      }
      if (cl1 == 'lmfd') {
        b = polm(-unclass(e1$b))
        return(lmfd(a = e1$a, b = b))
      }
      if (cl1 == 'rmfd') {
        d = polm(-unclass(e1$d))
        return(rmfd(c = e1$c, d = d))
      }
      if (cl1 == 'stsp') {
        return(stsp(A = e1$A, B = e1$B, C = -e1$C, D = -e1$D))
      }
      if (cl1 == 'pseries') {
        e1 = -unclass(e1)
        e1 = structure(e1, class = c('pseries', 'ratm'))
        return(e1)
      }
      if (cl1 == 'zvalues') {
        z = attr(e1,'z')
        e1 = -unclass(e1)
        e1 = structure(e1, z = z, class = c('zvalues', 'ratm'))
        return(e1)
      }
      stop('unsupported class: "',cl1,'"')
    }
    stop('unsupported unary operator: "',.Generic,'"')
  } 

  # power operator e1^n ###########################################
  if (.Generic == '^') {
    d1 = unname(dim(e1))
    cl1 = class(e1)[1]
    a2 = unclass(e2)
    
    if ( ( length(a2) != 1 ) || (a2 %% 1 != 0 ) ) {
      stop('unsupported power!')
    }
    
    if ( d1[1] != d1[2] ) {
      stop('power operation is only defined for non empty, square rational matrices!')
    }
    
    # a2 == 1 ######################
    if (a2 == 1) {
      return(e1)
    } # a2 = 1 
    
    # a2 == 0 ######################
    if (a2 == 0) {
      if (cl1 == 'polm') {
        return(polm(diag(d1[1])))
      }
      if (cl1 == 'lmfd') {
        return(lmfd(a = diag(d1[1]), b = diag(d1[1])))
      }
      if (cl1 == 'rmfd') {
        return(rmfd(c = diag(d1[1]), d = diag(d1[1])))
      }
      if (cl1 == 'stsp') {
        return(stsp(A = matrix(0, nrow = 0, ncol = 0), 
                    B = matrix(0, nrow = 0, ncol = d1[1]),
                    C = matrix(0, nrow = d1[1], ncol = 0),
                    D = diag(d1[1])))
      }
      if (cl1 == 'pseries') {
        if (d1[1] > 0) {
          e1 = unclass(e1)
          e1[,,-1] = 0
          e1[,,1] = diag(d1[1])
          e1 = structure(e1, class = c('pseries', 'ratm'))
        }
        return(e1)
      }
      if (cl1 == 'zvalues') {
        if (d1[1] > 0) {
          z = attr(e1,'z')
          e1 = array(diag(d1[1]), dim = d1)
          e1 = structure(e1, z = z, class = c('zvalues', 'ratm'))
        }
        return(e1)
      }
      stop('unsupported class: "',cl1,'"')
    } # a2 = 0
    
    # upgrade "lmfd" to "stsp" objects
    if (cl1 == 'lmfd') {
      e1 = as.stsp(e1)
      cl1 = 'stsp'
      d1 = unname(dim(e1))
    }
    
    # upgrade "rmfd" to "stsp" objects
    if (cl1 == 'rmfd') {
      e1 = as.stsp(e1)
      cl1 = 'stsp'
      d1 = unname(dim(e1))
    }
    
    # a2 > 1 ######################
    if (a2 > 1) {
      e = e1
      for (i in (2:a2)) {
        e = e %r% e1
      }
      return(e)
    } # a2 > 1
    
    # a2 < 0 ######################
    
    if (d1[1] <= 0) {
      stop('power operation with negative power is only defined for non empty, square rational matrices!')
    }
    
    # convert "polm" to "stsp" objects 
    if (cl1 == 'polm') {
      e1 = as.stsp(e1)
      cl1 = 'stsp'
      d1 = unname(dim(e1))
    }
    
    if (cl1 == 'stsp')  {
      # compute inverse 
      A = e1$A
      B = e1$B
      C = e1$C
      D = e1$D
      D = try(solve(D), silent = TRUE)
      if (inherits(D, 'try-error')) {
        stop('could not compute state space representation of inverse (D is singular)')
      }
      B = B %*% D
      e1 = stsp(A = A - B %*% C, B = B, C = -D %*% C, D = D)
      
      e = e1
      for (i in iseq(2,abs(a2))) {
        e = e %r% e1
      }
      return(e)
    }
    
    if (cl1 == 'pseries')  {
      # compute inverse 
      a = unclass(e1)
      m = dim(a)[1]           # we could also use d1!
      lag.max = dim(a)[3] - 1
      if (lag.max < 0) {
        # this should not happen?!
        stop('impulse response contains no lags!')
      }
      
      # b => inverse impulse response 
      b = array(0, dim = c(m,m,lag.max+1))
      
      # a[0] * b[0] = I
      b0 = try(solve(matrix(a[,,1], nrow = m, ncol = m)))
      if (inherits(b0, 'try-error')) {
        stop('impulse response is not invertible (lag zero coefficient is singular)')
      }
      b[,,1] = b0
      for (i in iseq(1,lag.max)) {
        # a[i] * b[0] + ... + a[0] b[i] = 0
        for (j in (1:i)) {
          b[,,i+1] = b[,,i+1] + matrix(a[,,j + 1], nrow = m, ncol = m) %*% 
                                matrix(b[,,i - j + 1], nrow = m, ncol = m)
        }
        b[,,i+1] = -b0 %*% matrix(b[,,i+1], nrow = m, ncol = m)
      }      
      
      # convert to pseries object
      class(b) = c('pseries','ratm')
      
      e = b
      for (i in iseq(2,abs(a2))) {
        e = e %r% b
      }
      return(e)
    }
    
    if (cl1 == 'zvalues')  {
      z = attr(e1, 'z')
      e1 = unclass(e1)
      # compute inverse 
      for (i in (1:length(z))) {
        ifr = try(solve(matrix(e1[,,i], nrow = d1[1], ncol = d1[1])), silent = TRUE)
        if (inherits(ifr, 'try-error')) {
          ifr = matrix(NA_real_, nrow = d1[1], ncol=d1[1])
        }
        e1[,,i] = ifr
      }
      e1 = structure(e1, z = z, class = c('zvalues', 'ratm'))
      
      e = e1
      for (i in iseq(2,abs(a2))) {
        e = e %r% e1
      }
      return(e)
    }
    
    # this should not happen!
    stop('unsupported class: "',cl1,'"')
  } 
  
  
  # elementwise operations '*', '+', '-', '%/%', '%%' ################################
  
  # make sure that both arguments have the same class!
  out = upgrade_objects(force = TRUE, e1, e2)
  e1 = out[[1]]
  e2 = out[[2]]

  d1 = unname(dim(e1))
  cl1 = class(e1)[1]
  gr1 = match(cl1, c('polm','lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)
  d2 = unname(dim(e2))
  cl2 = class(e2)[1]
  gr2 = match(cl1, c('polm','lmfd','rmfd','stsp','pseries','zvalues'), nomatch = 0)
  
  if (d1[1]*d1[2] == 1) {
    # e1 is a scalar 
    e1 = inflate_object(e1, m = d2[1], n = d2[2])
    d1 = unname(dim(e1))
  }
  if (d2[1]*d2[2] == 1) {
    # e2 is a scalar 
    e2 = inflate_object(e2, m = d1[1], n = d1[2])
    d2 = unname(dim(e2))
  }

  if (any(d1[1:2] != d2[1:2]))  {
    stop('the rational matrices e1, e2 are not compatible (nrow(e1) != nrow(e2) or ncol(e1) != ncol(e2))')
  }
  
  # __elementwise multiplication '*' ################################
  if (.Generic == '*') {
    # elementwise addition/substraction
    
    if (cl1 == 'polm') {
      e = polm(emult_poly(unclass(e1), unclass(e2)))
      e = prune(e, tol = 0)
      return(e)
    }
    
    if (cl1 == 'stsp') {
      e = emult_stsp(e1, e2)
      return(e)
    }
    
    if (cl1 == 'pseries') {
      e = emult_pseries(unclass(e1), unclass(e2))
      class(e) = c('pseries','ratm')
      return(e)
    }
    
    if (cl1 == 'zvalues') {
      z1 = attr(e1,'z')
      z2 = attr(e2,'z')
      if (!isTRUE(all.equal(z1, z2))) {
        stop('the complex points z of the two frequency responses do not coincide')
      }
      e1 = unclass(e1)
      e2 = unclass(e2)
      e = array(0, dim = c(d1[1], d2[2], length(z1)))
      for (i in (1:length(z1))) {
        e[,,i] = matrix(e1[,,i], nrow = d1[1], ncol = d1[2]) * matrix(e2[,,i], nrow = d2[1], ncol = d2[2])
      }
      
      e = structure(e, z = z1, class = c('zvalues', 'ratm'))
      return(e)
    }      
    
  } # elementwise multiplication '*' 
  
  # __elementwise addition/substraction '+', '-' ################################
  if ((.Generic == '+') || (.Generic == '-')) {
    # elementwise addition/substraction
    
    if (.Generic == '-') e2 = -e2

    if (cl1 == 'polm') {
      # polynomial matrices
      e1 = unclass(e1)
      e2 = unclass(e2)
      e = array(0, dim = c(d1[1], d1[2], max(d1[3], d2[3])+1))
      if (d1[3]>=0) e[,,1:(d1[3]+1)] = e1
      if (d2[3]>=0) e[,,1:(d2[3]+1)] = e[,,1:(d2[3]+1),drop = FALSE] + e2
      return(polm(e))
    }
    
    if (cl1 == 'stsp') {
      # statespace realization
      e1 = unclass(e1)
      e2 = unclass(e2)
      
      A = bdiag(e1[iseq(1,       d1[3]),       iseq(1,       d1[3])        , drop = FALSE],             
                e2[iseq(1,       d2[3]),       iseq(1,       d2[3])        , drop = FALSE])
      B = rbind(e1[iseq(1,       d1[3]),       iseq(d1[3]+1, d1[3]+d1[2])  , drop = FALSE], 
                e2[iseq(1,       d2[3]),       iseq(d2[3]+1, d2[3]+d2[2])  , drop = FALSE])
      C = cbind(e1[iseq(d1[3]+1, d1[3]+d1[1]), iseq(1,       d1[3])        , drop = FALSE], 
                e2[iseq(d2[3]+1, d2[3]+d2[1]), iseq(1,       d2[3])        , drop = FALSE])
      D =       e1[iseq(d1[3]+1, d1[3]+d1[1]), iseq(d1[3]+1, d1[3]+d1[2])  , drop = FALSE]     + 
        e2[iseq(d2[3]+1, d2[3]+d2[1]), iseq(d2[3]+1, d2[3]+d2[2])  , drop = FALSE]
      return(stsp(A = A, B = B, C = C, D = D))
    }
    
    if (cl1 == 'pseries') {
      # print(attributes(e1))
      # print(attributes(e2))
      # impulse response
      e1 = unclass(e1)
      e2 = unclass(e2)
      if (d1[3] <= d2[3]) {
        e1 = e1 + e2[,,iseq(1, d1[3]+1), drop = FALSE]
        e1 = structure(e1, class = c('pseries', 'ratm'))
        return(e1)
      }
      e2 = e2 + e1[,,iseq(1, d2[3]+1), drop = FALSE]
      e2 = structure(e2, class = c('pseries', 'ratm'))
      return(e2)
    }
    
    if (cl1 == 'zvalues') {
      # frequency response
      z1 = attr(e1,'z')
      z2 = attr(e2,'z')
      if (!isTRUE(all.equal(z1, z2))) {
        stop('the complex points z of the two frequency responses do not coincide')
      }
      e1 = unclass(e1)
      e2 = unclass(e2)
      e1 = e1 + e2
      e1 = structure(e1, z = z1, class = c('zvalues', 'ratm'))
      return(e1)
    }

  }  # elementwise addition/substraction '+', '-' 
  

  # __elementwise (polynomial) division #########################################
  if (.Generic == '%/%') {
    if ( cl1 != 'polm' ) {
      stop('elementwise (polynomial) divsision "%/%" is only implemented for "polm" objects')
    }
    
    e = polm(array(0, dim = c(d1[1], d1[2], 1)))
    if (d1[1]*d1[2] == 0) {
      return( e )
    }
    
    a1 = unclass(e1)
    a2 = unclass(e2)
    for ( i in (1:d1[1]) ) {
      for (j in (1:d1[2])) {
        #        print(str(a))
        #        print(polm_div(a1[i,j,], a2[i,j,]))
        e[i,j] = poly_div(a1[i,j,], a2[i,j,])
      }
    }
    # skip leading zero coefficient matrices
    e = prune(e, tol = 0)
    return( e )
  }
  
  # __elementwise (polynomial) remainder #########################################
  if (.Generic == '%%') {
    if ( cl1 != 'polm' ) {
      stop('elementwise (polynomial) remainder "%%" is only implemented for "polm" objects')
    }
    
    e = polm(array(0, dim = c(d1[1], d1[2], 1)))
    
    if (d1[1]*d1[2] == 0) {
      return( e )
    }
    
    a1 = unclass(e1)
    a2 = unclass(e2)
    for ( i in (1:d1[1]) ) {
      for (j in (1:d1[2])) {
        #        print(str(a))
        #        print(polm_div(a1[i,j,], a2[i,j,]))
        e[i,j] = poly_rem(a1[i,j,], a2[i,j,])
      }
    }
    # skip leading zero coefficient matrices
    e = prune(e, tol = 0)
    return( e )
  }

  stop('unsupported operator: "',.Generic,'"')
}
