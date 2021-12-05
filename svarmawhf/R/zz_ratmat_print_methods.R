# print.____ methods ##############################################################

#' Print Methods
#'
#' These functions were originally part of the R-package \strong{rationalmatrices}.
#' 
#' @param x rational matrix object, i.e. a \code{\link{polm}}, \code{\link{lmfd}}, 
#'          \code{\link{stsp}}, \code{\link{pseries}} or \code{\link{zvalues}} object.
#' @param digits (integer) if non \code{NULL} then correspondingly rounded numbers are printed, 
#'        see \code{\link{round}}.
#' @param format (character string) selects specific output formats. Note that 
#'        \code{\link{stsp}} objects have no format option. The option \code{'character'} 
#'        is only implemented for polynomials and LMFDs with real coefficients.
#' @param ... Further parameters are ignored.
#'
#' @return \code{invisible(x)}
#'
#' @rdname print
#' @name print methods
#' 
#' @examples
#' # for polynomials six different print formats are implemented ###################
#' a = test_polm(dim = c(2,3), degree = 2, random = TRUE)
#' 
#' for (fmt in c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character")) {
#'    cat('\nformat =', fmt, '\n')
#'    print(a, digits = 2, format = fmt)
#' }
#' 
#' # "empty" (2 x 0) polynomial matrix (degree = 2)
#' a = test_polm(dim = c(2,0), degree = 0)
#' print(a)
#' 
#' # random (2 x 1) polynomial matrix with complex coefficients (degree = 2)
#' a = polm(array(complex(real = stats::rnorm(2*1*3), 
#'                        imaginary = stats::rnorm(2*1*3)), dim = c(2,1,3)))
#' print(a, digits = 2)
#' \dontrun{
#' # the format option 'character' is only implemented for polynomials matrices 
#' # with real coefficients!
#' print(a, digits = 2, format = 'character')
#' }
#'
#' # print a rational matrix in 'lmfd' form 
#' a = test_lmfd(dim = c(2,3), degrees = c(2,1))
#' print(a, digits = 2, format = 'character')
#' 
NULL


# internal function 
# convert univariate polynomial to string
# c vector of doubles. (No complex numbers allowed)
as_character_poly = function(c) {
  if (is.complex(c)) {
    stop('only the real case is implemented')
  }
  if (length(c) > 0) {
    c = c[rev(cumprod(rev(c == 0))) == 0]
  }
  if (length(c) == 0) {
    return('0')
  }
  
  p = length(c)-1
  powers = (0:p)

  # skip zero coefficients   
  non_zero = (c != 0)
  c = c[non_zero]
  powers = powers[non_zero]
  
  powers = paste('z^', powers, sep='')
  powers[powers == 'z^0'] = ''
  powers[powers == 'z^1'] = 'z'
  
  signs = ifelse(c < 0, '- ', '+ ')
  signs[1] = ifelse(signs[1] == '- ', '-', '')
  
  c = as.character(abs(c))
  c[ (c == '1') & (powers != '') ] = ''
  
  paste(signs, c, powers, sep='', collapse = ' ')
}

# internal function
# print 3D array
print_3D = function(a, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character')) {
  dim = dim(a)
  m = dim[1]
  n = dim[2]
  p = dim[3]
  # empty array -> do nothing
  if (min(dim) == 0) return(invisible(NULL))
  
  # a must have full 'dimnames'
  names = dimnames(a)
  inames = names[[1]]
  jnames = names[[2]]
  znames = names[[3]]
  
  # round 
  if (!is.null(digits)) a = round(a, digits)
  
  format = match.arg(format)
  
  if (format == 'character') {

    # convert vector of coefficients to character representation of a polynomial
    a = apply(a, MARGIN = c(1,2), FUN = as_character_poly)
    
    # add column names (jnames)
    a = rbind(jnames, a)
    
    # add row names (inames)
    a = cbind( c('',inames), a)
    
    # right justify columns 
    w = nchar(a)
    w = apply(w, MARGIN = 2, FUN = max)
    for (j in (1:(n+1))) {
      fmt = paste('%', w[j], 's', sep='')
      pad = function(s) { sprintf(fmt, s) }
      a[,j] = apply(a[,j,drop = FALSE], MARGIN = 1, FUN = pad)
    }

    # convert matrix a to a string
    a = apply(a, MARGIN = 1, FUN = paste, collapse = '  ')
    a = paste(a, collapse = '\n')
    cat(a,'\n')
  }
  
  
  if (format == 'i|jz') {
      # create a vector of the form 
      # j[1],...,j[n],j[1],...,j[n],...
      jnames = rep(jnames, p)
      # create a vector of the form 
      # z[1],'',...,'',z[2],'',...,'',
      if (n > 1) {
        znames = as.vector(rbind(znames, 
                                 matrix('', nrow = n-1, ncol = p)))
      }
      
      dim(a) = c(m,n*p)
      rownames(a) = inames
      colnames(a) = paste(znames, jnames, sep = ' ')
      print(a)  
  }
  
  if (format == 'i|zj') {
    # create a vector of the form 
    # z[1],...,z[p],z[1],...,z[p],...
    znames = rep(znames, n)
    # create a vector of the form 
    # j[1],'',...,'',j[2],'',...,'',
    if (p > 1) {
      jnames = as.vector(rbind(jnames, 
                               matrix('', nrow = p-1, ncol = n)))
    }
    
    a = aperm(a, c(1,3,2))
    dim(a) = c(m,p*n)
    rownames(a) = inames
    colnames(a) = paste(jnames, znames, sep = ' ')
    print(a)  
  }
  
  if (format == 'iz|j')  {
    # create a vector of the form 
    # i[1],...,i[m],i[1],...,i[m],...
    inames = rep(inames, p)
    # create a vector of the form 
    # z[1],'  ',...,'  ',z[2],'  ',...,'  ',
    if (m > 1) {
      znames = as.vector(rbind( znames, 
                                matrix(' ', nrow = m-1, ncol = p)))
    }
    # right justify
    fmt = paste('%', max(nchar(znames)), 's', sep='')
    pad = function(s) { sprintf(fmt, s) }
    znames = as.vector(apply(matrix(znames, ncol = 1), MARGIN = 1, FUN = pad))
    
    a = aperm(a, c(1,3,2))
    dim(a) = c(m*p, n)
    rownames(a) = paste(znames, inames, sep = ' ')
    colnames(a) = jnames
    print(a)
  }
  
  if (format == 'zi|j')  {
    # create a vector of the form 
    # z[1],...,z[p],z[1],...,z[p],...
    znames = rep(znames, m)
    # create a vector of the form 
    # i[1],'  ',...,'  ',i[2],'  ',...,'  ',
    if (p > 1) {
      inames = as.vector(rbind( inames, 
                                matrix(' ',nrow = p-1, ncol = m)))
    }
    # right justify
    fmt = paste('%', max(nchar(inames)), 's', sep='')
    pad = function(s) { sprintf(fmt, s) }
    inames = as.vector(apply(matrix(inames, ncol = 1), MARGIN = 1, FUN = pad))

    a = aperm(a, c(3,1,2))
    dim(a) = c(p*m, n)
    rownames(a) = paste(inames, znames, sep = ' ')
    colnames(a) = jnames
    print(a)
  }
  
  if (format == 'i|j|z') {
    # the last case 'i|j|z' just uses the R default print of 3D array
    print(a)
  }
    
  return(invisible(NULL))
}


#' @rdname print
#' @export
print.polm = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)
  
  a = unclass(x)
  m = dim(a)[1]
  n = dim(a)[2]
  p = dim(a)[3]-1
  
  cat('(',m,'x',n,') matrix polynomial with degree <=', p,'\n')
  if ((m*n*(p+1)) == 0) {
    return(invisible(x))
  }
  
  # if (format == 'character') {
  #   # rounding digits
  #   if (!is.null(digits)) {
  #     a = round(a, digits)
  #   }
  #   
  #   a = apply(a, MARGIN = c(1,2), FUN = as_character_poly)
  #   a = rbind(paste('[,',1:n,']', sep = ''), a)
  #   w = nchar(a)
  #   w = apply(w, MARGIN = 2, FUN = max)
  #   for (j in (1:n)) {
  #     fmt = paste('%',w[j],'s',sep='')
  #     pad = function(s) { sprintf(fmt, s) }
  #     a[,j] = apply(a[,j,drop = FALSE], MARGIN = 1, FUN = pad)
  #   }
  #   rnames = matrix(c('',paste('[',1:m,',]', sep = '')), ncol = 1)
  #   w = max(nchar(rnames))
  #   fmt = paste('%', w, 's', sep='')
  #   pad = function(s) { sprintf(fmt, s) }
  #   rnames = apply(rnames, MARGIN = 1, FUN = pad)
  #   a = cbind(rnames, a)
  #   a = apply(a, MARGIN = 1, FUN = paste, collapse = ' ')
  #   a = paste(a, collapse = '\n')
  #   cat(a,'\n')
  # }
  
  if ((format == 'character') && (is.complex(a))) {
    stop('the format option "character" is only implemented for polynomials with real coefficients')
  }
  
  # use the above defined internal function print_3D
  dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                     paste('[,', 1:n, ']', sep = ''),
                     paste('z^',0:p, sep = ''))
  print_3D(a, digits, format)
  
  invisible(x)
}


#' @rdname print
#' @export
print.lmfd = function(x, digits = NULL, format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z','character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)
  
  d = attr(x, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]
  
  cat('( ', m, ' x ', n,' ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = ', 
      p, ', q = ', q, ')\n', sep = '')
  
  if ((format == 'character') && (is.complex(unclass(x)))) {
    stop('the format option "character" is only implemented for LMFDs with real coefficients')
  }

  if ((m*m*(p+1)) > 0) {
    cat('left factor a(z):\n')
    
    a = unclass(x$a)
    
    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:m, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(a, digits, format)
  }
  
  if ((m*n*(q+1)) > 0) {
    cat('right factor b(z):\n')
    
    a = unclass(x$b)
    
    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(a, digits, format)
  }
  
  invisible(x)
}
