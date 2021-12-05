# print.____ methods ##############################################################

#' Coerce Scalar Polynomials to Character Strings
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' \cr
#' This utility coerces a scalar polynomial (given by the vector of coefficients)
#' to a character string. The following "formats" are implemented.
#' \code{syntax = "txt"} returns a simple text representation,
#' \code{syntax = "TeX"} renders the coefficients to string in "TeX" syntax and
#' \code{syntax = "expression"} gives a string which may be rendered to
#' an \code{R} expression with \code{\link[base]{parse}}. This expression
#' may be used to evaluate the polynomial and for annotating plots,
#' see \code{\link[grDevices]{plotmath}} and the examples below.
#'
#' @param coefs (numeric) vector of coefficients.
#' @param syntax (character string) determines the format of the output string.
#' @param x (character string) polynomial variable.
#'
#' @return character string.
#' @export
#'
#' @examples
#' coefs = c(1, 2.3, 0, -1, 0)
#'
#' as_txt_scalarpoly(coefs, syntax = 'txt', x = 'x')
#' as_txt_scalarpoly(coefs, syntax = 'TeX', x = '\\alpha')
#' as_txt_scalarpoly(coefs, syntax = 'expression', x = 'z')
#'
#' \dontrun{
#' # the case syntax = "expression" may be used e.g. as follows
#'
#' # make_polyfun creates a "closure" which evaluates the polynomial at given points
#' # note that this simple version does not work for zero polnomials!
#' make_polyfun = function(coefs) {
#'   expr = parse(text = as_txt_scalarpoly(coefs, 'expression', 'x'))
#'   fun = function(x) {
#'     return(eval(expr))
#'   }
#'   return(fun)
#' }
#'
#' a = make_polyfun(coefs)
#' a(1)   # return the value  of the polynomial at x = 1
#' a(1:5) # return the values of the polynomial at x = 1,2,3,4,5
#'
#' # create a plot
#' x_grid = seq(from = -1, to = 1, length.out = 101)
#' plot(x_grid, a(x_grid), type = 'l', xlab = 'x', ylab = 'a(x)',
#'      main = parse(text = paste('a(x) == ',
#'           as_txt_scalarpoly(coefs, syntax = 'expression', x = 'x'))))
#' }
as_txt_scalarpoly = function(coefs, syntax = c('txt', 'TeX', 'expression'),
                                   x = 'z') {
  coefs = as.vector(coefs)
  if (!is.numeric(coefs)) {
    stop('"coefs" must be a numeric vector')
  }
  syntax = match.arg(syntax)

  if ((length(coefs) == 0) || all(coefs == 0)) {
    return('0')
  }
  powers = (0:(length(coefs)-1))

  # skip zero coefficients
  non_zero = (coefs != 0)
  coefs = coefs[non_zero]
  powers = powers[non_zero]

  # convert powers to character strings
  if (syntax == 'txt') {
    # x^k
    powers_txt = paste(x, '^', powers, sep = '')
  } else {
    # x^{k}
    powers_txt = paste(x, '^{', powers, '}', sep = '')
    fmt = 'x^{k}'
  }
  powers_txt[powers == 0] = ''
  powers_txt[powers == 1] = x
  powers = powers_txt

  signs = ifelse(coefs < 0, '- ', '+ ')
  signs[1] = ifelse(coefs[1] < 0, '-', '')

  # convert coefficients to character strings
  coefs = paste(abs(coefs))
  coefs[ (coefs == '1') & (powers != '') ] = ''

  if (syntax == 'expression') {
    mults = rep('*', length(coefs))
    mults[ (coefs == '') | (powers == '') ] = ''
  } else {
    mults = rep('', length(coefs))
  }

  txt = paste(signs, coefs, mults, powers, sep = '', collapse = ' ')
  return(txt)

}

#' Coerce Scalar Polynomial Filters to Character Strings
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' \cr
#' This utility coerces a scalar polynomial filter (given by the vector
#' of coefficients) to a character string. The following "formats" are implemented.
#' \code{syntax = "txt"} returns a simple text representation,
#' \code{syntax = "TeX"} renders the coefficients to string in "TeX" syntax and
#' \code{syntax = "expression"} gives a string which may be rendered to
#' an \code{R} expression with \code{\link[base]{parse}}. This expression
#' may be used to evaluate the filter and for annotating plots,
#' see \code{\link[grDevices]{plotmath}} and the examples below.
#'
#' @param coefs (numeric) vector of coefficients.
#' @param syntax (character string) determines the format of the output string.
#' @param x (character string) names the "input" series.
#' @param t (character string) names the "time-index".
#'
#' @return character string.
#' @export
#'
#' @seealso \code{\link{as_txt_scalarpoly}} and \code{\link{as_tex_matrixfilter}}.
#'
#' @examples
#' coefs = c(1, 2.3, 0, -1, 0)
#'
#' as_txt_scalarfilter(coefs, syntax = 'txt', x = 'x', t = 't')
#' as_txt_scalarfilter(coefs, syntax = 'TeX', x = 'x', t = 's')
#' as_txt_scalarfilter(coefs, syntax = 'expression', x = 'x', t = 'k')
#'
#' \dontrun{
#' # the case syntax = "expression" may be used e.g. as follows
#'
#' # make_filterfun creates a "closure" which computes the filter-output
#' # note that this simple version does not work for zero filters!
#' make_filterfun = function(coefs) {
#'   p = length(coefs) - 1
#'   expr = parse(text = as_txt_scalarfilter(coefs, 'expression', 'x', 't'))
#'   fun = function(x, t) {
#'     # x, t must be vectors
#'     y = rep(NA_real_, length(t))
#'     t0 = t
#'     y = rep(NA_real_, length(t))
#'
#'     i = ((t0 > p) & (t0 <= length(x)))
#'     t = t0[i]
#'     if (any(i)) y[i] = eval(expr)
#'
#'     return(y)
#'   }
#'   return(fun)
#' }
#'
#' coefs = rep(1, 4) / 4  # represents a moving average of length 4.
#' a = make_filterfun(coefs)
#' u = rnorm(100)       # input series
#' a(u, 1)    # return the value of the output series at t = 1
#'            # this value is not defined due to missing initial values
#' a(u, 1:10) # return the values of the output series at t = 1,..,10
#'
#' # create a plot
#' plot(1:length(u), u, type = 'n', xlab = 'time', ylab = '')
#' grid()
#' lines(1:length(u), u, col = 'black', lwd = 1)
#' lines(1:length(u), a(u, 1:length(u)), col = 'red', lwd = 2)
#' legend('topright', bty = 'n',
#'        fill = c('black', 'red'),
#'        legend = c(expression(u[t]),
#'                   parse(text = paste('x[t] == ',
#'                      as_txt_scalarfilter(coefs, 'expression', 'u','t')))) )
#' }
as_txt_scalarfilter = function(coefs, syntax = c('txt', 'TeX', 'expression'),
                               x = 'z', t = 't') {
  coefs = as.vector(coefs)
  if (!is.numeric(coefs)) {
    stop('"coefs" must be a numeric vector')
  }
  syntax = match.arg(syntax)

  if ((length(coefs) == 0) || all(coefs == 0)) {
    return('0')
  }
  lags = (0:(length(coefs)-1))

  # skip zero coefficients
  non_zero = (coefs != 0)
  coefs = coefs[non_zero]
  lags = lags[non_zero]

  # convert lags to character strings
  if (syntax == 'TeX') {
    # x_{t-k}
    lags_txt = paste(x, '_{', t, '-', lags, '}', sep = '')
    lags_txt[lags == 0] = paste(x, '_{', t, '}', sep = '')
  } else {
    # x[t-k]
    lags_txt = paste(x, '[', t, '-', lags, ']', sep = '')
    lags_txt[lags == 0] = paste(x, '[', t, ']', sep = '')
  }
  lags = lags_txt

  signs = ifelse(coefs < 0, '- ', '+ ')
  signs[1] = ifelse(coefs[1] < 0, '-', '')

  # convert coefficients to character strings
  coefs = paste(abs(coefs))
  coefs[ (coefs == '1') & (lags != '') ] = ''

  if (syntax == 'expression') {
    mults = rep('*', length(coefs))
    mults[ (coefs == '') | (lags == '') ] = ''
  } else {
    mults = rep('', length(coefs))
  }
  txt = paste(signs, coefs, mults, lags, sep = '', collapse = ' ')
  return(txt)

}

#' TeX Matrix
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' 
#' @param x matrix, where \code{paste(x[i,j])} returns a valid "TeX" string.
#'
#' @return character string.
#' @export
#'
#' @examples
#' as_tex_matrix(diag(1:2, nrow = 2, ncol = 3))
as_tex_matrix = function(x) {
  if ( !is.matrix(x) ) stop('"x" must be a matrix')

  m = nrow(x)
  n = ncol(x)

  if (length(x) == 0) return('\\begin{pmatrix}\n\\end{pmatrix}')

  tex = '\\begin{pmatrix}\n'
  for (i in (1:m)) {
    tex = paste(tex, paste(x[i,], collapse = ' & '), '\\\\\n', sep = '  ')
  }
  tex = paste(tex, '\\end{pmatrix}', sep = '')
  return(tex)
}


#' TeX Matrix Polynomials
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' 
#' @param coefs 3-dimensional array with the coefficients of the matrix polynomial.
#' @param x (character string) polynomial variable or process variable.
#' @param as_matrix_of_polynomials boolean.
#'
#' @return character string.
#' @export
#'
#' @examples
#' coefs = array(round(rnorm(2*3*1), 1), dim = c(2,3,2))
#'
#' as_tex_matrixpoly(coefs)
#' as_tex_matrixpoly(coefs, x = 'x', as_matrix_of_polynomials = FALSE)
as_tex_matrixpoly = function(coefs, x = 'z', as_matrix_of_polynomials = TRUE) {
  # only some basic checks
  if ( (!is.array(coefs)) || (length(dim(coefs)) != 3) || (!is.numeric(coefs)) ) {
    stop('"coefs" must be 3-dimensional numeric array')
  }

  d = dim(coefs)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  if ((m*n) == 0) {
    return('\\begin{pmatrix}\n\\end{pmatrix}')
  }

  if ((m*n) == 1) {
    return(as_txt_scalarpoly(coefs, syntax = 'TeX', x = x))
  }

  if ((p < 0) || all(coefs == 0)) {
    return(as_tex_matrix(matrix(0, nrow = m, ncol = n)))
  }

  if (as_matrix_of_polynomials) {
    tex = apply(coefs, MARGIN = c(1,2), FUN = as_txt_scalarpoly,
                syntax = 'TeX', x = x)
    tex = as_tex_matrix(tex)
    return(tex)
  }

  # print as polynomial with matrix coefficients

  powers = (0:p)
  # coerce powers to character strings of the form x^{k}
  powers_txt = paste(x, '^{', powers, '}', sep = '')
  powers_txt[powers == 0] = ''
  powers_txt[powers == 1] = x
  powers = powers_txt

  tex = ''
  for (k in (0:p)) {
    a = matrix(coefs[,,k+1], nrow = m, ncol = n)
    if ( !all(a == matrix(0, nrow = m, ncol = n)) ) {
      # non-zero coefficient matrix
      if (tex != '' ) tex = paste(tex, '+\n')
      if ( (m == n) && all(a == diag(m)) ) {
        # coefficient matrix is identity matrix
          tex = paste(tex, ' I_{', m, '} ', powers[k+1], sep = '')
      } else {
        tex = paste(tex, as_tex_matrix(a), powers[k+1])
      }
    }
  }

  return(tex)
}


#' TeX Matrix Polynomial Filters
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' 
#' @param coefs 3-dimensional array with the coefficients of the filter.
#' @param x (character string) polynomial variable or process variable.
#' @param t (character string) time/index variable.
#'
#' @return character string.
#' @export
#'
#' @examples
#' coefs = array(round(rnorm(2*3*1), 1), dim = c(2,3,2))
#'
#' as_tex_matrixfilter(coefs, x = '\\epsilon', t = 's')
as_tex_matrixfilter = function(coefs, x = 'z', t = 't') {
  # only some basic checks
  if ( (!is.array(coefs)) || (length(dim(coefs)) != 3) || (!is.numeric(coefs)) ) {
    stop('"coefs" must be 3-dimensional numeric array')
  }

  d = dim(coefs)
  m = d[1]
  n = d[2]
  p = d[3] - 1

  if ((m*n) == 0) {
    return('\\begin{pmatrix}\n\\end{pmatrix}')
  }

  if ((m*n) == 1) {
    return(as_txt_scalarfilter(coefs, syntax = 'TeX', x = x, t = t))
  }

  if ((p < 0) || all(coefs == 0)) {
    tex = as_tex_matrix(matrix(0, nrow = m, ncol = n))
    return( paste(tex, ' ', x, '_{', t, '}', sep = '') )
  }

  lags = (0:p)
  # coerce lags to character strings of the form x_{t-k}
  lags_txt = paste(x, '_{', t, '-', lags, '}', sep = '')
  lags_txt[lags == 0] = paste(x, '_{', t, '}', sep = '')
  lags = lags_txt

  tex = ''
  for (k in (0:p)) {
    a = matrix(coefs[,,k+1], nrow = m, ncol = n)
    if ( !all(a == matrix(0, nrow = m, ncol = n)) ) {
      # non-zero coeffient matrix
      if (tex != '' ) tex = paste(tex, '+\n')
      if ( (m==n) && all(a == diag(m)) ) {
        # coefficient matrix is identity matrix
        tex = paste(tex, ' I_{', m, '} ', lags[k+1], sep = '')
      } else {
        tex = paste(tex, as_tex_matrix(a), lags[k+1])
      }
    }
  }

  return(tex)
}

#' Print 3D array
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' 
#' @param a 3-dimensional (numeric or complex) array
#' @param digits (integer) if non \code{NULL} then correspondingly rounded numbers
#'        are printed, see \code{\link{round}}.
#' @param format format string
#'
#' @return invisible(NULL)
#' @export
#' @keywords internal
print_3D = function(a, digits = NULL,
                    format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z', 'character')) {
  # no checks!!!!

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

  if (format == 'character') {

    # convert vector of coefficients to character representation of a polynomial
    a = apply(a, MARGIN = c(1,2), FUN = as_txt_scalarpoly, syntax = 'txt', x = 'z')

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

  return(invisible(NULL))
}

#' Print Methods
#'
#' This function was originally part of the R-package \strong{RLDM}.
#' 
#' @param x RLDM object, i.e. a \code{\link{armamod}}, \code{\link{stspmod}},
#'          \code{\link{impresp}}, \code{\link{autocov}}, \code{\link{freqresp}},
#'          \code{\link{spectrum}} or \code{\link{fevardec}} object.
#' @param digits (integer) if non \code{NULL} then correspondingly rounded numbers are printed,
#'        see \code{\link{round}}.
#' @param format (character string) selects specific output formats. Note that
#'        \code{\link{stsp}} and \code{\link{fevardec}} objects have no format option.
#'        The option \code{'character'} is only implemented for (V)ARMA models.
#' @param ... Further parameters are ignored.
#'
#' @return \code{invisible(x)}
#'
#' @rdname print
#' @name print methods
#'
#' @examples
#' # for VARMA models six different print formats are implemented ###################
#' m = armamod(test_lmfd(dim = c(2,2), degrees = c(1,1)), sigma_L = diag(2))
#' print(m, digits = 2, format = "i|jz")
NULL

#' @rdname print
#' @export
print.armamod = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z', 'character'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  format = match.arg(format)
  if ((format == 'character') && (is.complex(unclass(x$sys)))) {
    stop('the format option "character" is only implemented for ARMA models with real coefficients.')
  }

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  d = attr(x$sys, 'order')
  m = d[1]
  n = d[2]
  p = d[3]
  q = d[4]

  cat(label, 'ARMA model [',d[1],',',d[2],'] with orders p = ', d[3], ' and q = ', d[4], '\n', sep = '')

  if ((m*m*(p+1)) > 0) {
    cat('AR polynomial a(z):\n')

    a = unclass(x$sys$a)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:m, ']', sep = ''),
                       paste('z^',0:p, sep = ''))
    print_3D(a, digits, format)
  }

  if ((m*n*(q+1)) > 0) {
    cat('MA polynomial b(z):\n')

    a = unclass(x$sys$b)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('z^',0:q, sep = ''))
    print_3D(a, digits, format)
  }

  if (n > 0) {
    cat('Left square root of noise covariance Sigma:\n')

    a = x$sigma_L
    dimnames(a) = list(paste('u[',1:n,']',sep = ''),paste('u[',1:n,']',sep = ''))
    print(a)
  }

  return(invisible(x))
}

#' @rdname print
#' @export
print.stspmod = function(x, digits = NULL, ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }

  d = attr(x$sys, 'order')
  m = d[1]
  n = d[2]
  s = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  cat(label, 'Statespace model [', m, ',', n, '] with s = ', s, ' states\n', sep = '')

  a = unclass(x$sys)
  attr(a, 'order') = NULL
  if (length(a) == 0) {
    return(invisible(x))
  }

  # rounding digits
  if (!is.null(digits)) {
    a = round(a, digits)
  }

  snames = character(s)
  if (s > 0) snames = paste('s[',1:s,']',sep = '')
  xnames = character(m)
  if (m > 0) {
    if ( !is.null(x$names) && is.character(x$names) && is.vector(x$names) && (length(x$names) == m) ) {
      xnames = x$names
    } else {
      xnames = paste('x[',1:m,']',sep = '')
    }
  }
  unames = character(n)
  if (n > 0) unames = paste('u[',1:n,']',sep = '')

  rownames(a) = c(snames, xnames)
  colnames(a) = c(snames, unames)
  print(a)

  if (n > 0) {
    cat('Left square root of noise covariance Sigma:\n')

    a = x$sigma_L
    dimnames(a) = list(paste('u[',1:n,']',sep = ''),paste('u[',1:n,']',sep = ''))
    print(a)
  }

  invisible(x)
}

#' @rdname print
#' @export
print.impresp = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$irf)
  m = d[1]
  n = d[2]
  lag.max = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  orth = FALSE
  if (n > 0) {
    sigma = x$sigma_L
    sigma = sigma %*% t(sigma)
    orth = isTRUE(all.equal(x$sigma_L, diag(n)))
  }

  if (orth) {
    cat(label, 'Orthogonalized impulse response [', m, ',', n, '] with ', lag.max, ' lags\n', sep = '')
  } else {
    cat(label, 'Impulse response [', m, ',', n, '] with ', lag.max, ' lags\n', sep = '')
  }

  if ((m*n*(lag.max+1)) > 0) {
    a = unclass(x$irf)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:n, ']', sep = ''),
                       paste('lag=', 0:lag.max, sep = ''))
    print_3D(a, digits, format)
  }

  invisible(x)
}

#' @rdname print
#' @export
print.autocov = function(x, digits = NULL,
                         format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$acf)
  m = d[1]
  lag.max = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  if (x$type == 'covariance') {
    label = paste(label, 'Autocovariance function [',d[1],',',d[2],'] with ', d[3], ' lags', sep = '')
  }
  if (x$type == 'correlation') {
    label = paste(label, 'Autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags', sep = '')
  }
  if (x$type == 'partial') {
    label = paste(label, 'Partial autocorrelation function [',d[1],',',d[2],'] with ', d[3], ' lags', sep = '')
  }

  n.obs = x$n.obs
  if (!is.null(n.obs)) {
    n.obs = as.integer(n.obs)[1]
  } else {
    n.obs = Inf
  }
  if ( is.finite(n.obs) ) {
    label = paste(label, ', sample size is ', n.obs, sep = '')
  }

  cat(label, '\n')

  if ((m*(lag.max+1)) > 0) {
    a = unclass(x$acf)

    # use the above defined internal function print_3D
    dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                       paste('[,', 1:m, ']', sep = ''),
                       paste('lag=', 0:lag.max, sep = ''))
    print_3D(a, digits, format)
  }
  invisible(x)
}




#' @rdname print
#' @export
print.fevardec = function(x, digits = NULL,
                          format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  vd = x$vd
  d = dim(vd)
  n = d[1]
  h.max = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')
  cat(label, 'Forecast error variance decompositon [', n, ',', n,'] maximum horizon = ', h.max, '\n', sep = '')

  if ((n*h.max)> 0) {
    names = x$names
    if (is.null(names)) {
      names = paste('y[', 1:n, ']', sep = '')
    }
    unames = paste('u[', 1:n, ']', sep = '')
    hnames = paste('h=', 1:h.max, sep='')

    dimnames(vd) = list(names, unames, hnames)
    print_3D(vd, digits, format)
  }

  invisible(x)
}


#' @rdname print
#' @export
print.freqresp = function(x, digits = NULL,
                          format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$frr)
  m = d[1]
  n = d[2]
  n.f = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  cat(label, 'Frequency response [',d[1],',',d[2],'] with ', d[3], ' frequencies\n', sep = '')

  if ((m*n*n.f) > 0) {
    a = unclass(x$frr)
    attr(a, 'z') = NULL
    f = (0:(n.f-1))/n.f

    # use the above defined internal function print_3D
    if ((format == 'i|jz') || (format == 'i|zj')) {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f[',1:n.f, ']', sep = ''))
    } else {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f=', round(f,3), sep = ''))
    }
    print_3D(a, digits, format)
  }

  invisible(x)
}


#' @rdname print
#' @export
print.spectrald = function(x, digits = NULL,
                           format = c('i|jz', 'i|zj', 'iz|j', 'zi|j', 'i|j|z'), ...) {
  if (!is.null(digits)) {
    digits = as.vector(as.integer(digits))[1]
  }
  format = match.arg(format)

  d = dim(x$spd)
  m = d[1]
  n = d[2]
  n.f = d[3]

  label = ''
  if (!is.null(x$label)) label = paste(x$label, ': ', sep = '')

  label = paste(label, 'Frequency response [',d[1],',',d[2],'] with ', d[3], ' frequencies', sep = '')

  n.obs = x$n.obs
  if (!is.null(n.obs)) {
    n.obs = as.integer(n.obs)[1]
  } else {
    n.obs = Inf
  }
  if ( is.finite(n.obs) ) {
    label = paste(label, ', sample size is ', n.obs, sep = '')
  }

  cat(label, '\n')

  if ((m*n.f) > 0) {
    a = unclass(x$spd)
    attr(a, 'z') = NULL
    f = (0:(n.f-1))/n.f

    # use the above defined internal function print_3D
    if ((format == 'i|jz') || (format == 'i|zj')) {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f[',1:n.f, ']', sep = ''))
    } else {
      dimnames(a) = list(paste('[', 1:m, ',]', sep = ''),
                         paste('[,', 1:n, ']', sep = ''),
                         paste('f=', round(f,3), sep = ''))
    }
    print_3D(a, digits, format)
  }

  invisible(x)
}

