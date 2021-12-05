# plot.____ methods ##############################################################


#' Plot Methods
#'
#' Plot impulse and frequency response functions.
#'
#' The parameter \code{which} determines what to plot in the case of \code{zvalues} objects.
#' \describe{
#' \item{modulus}{plot the moduli \code{abs(x[i,j])} versus frequencies \code{-Arg(z)/(2*pi)}.}
#' \item{phase}{plot the arguments \code{Arg(x[i,j])} versus frequencies \code{-Arg(z)/(2*pi)}.}
#' \item{nyquist}{plot the imaginary part \code{Im(x[i,j])} versus the real part \code{Re(x[i,j])}.}
#' \item{real}{plot the real part \code{Re(x[i,j])} versus the real part \code{Re(z)}. }
#' }
#' The choices \code{which = 'modulus'}, \code{which = 'phase'} and \code{which = 'nyquist'} 
#' are typically used for the case, where the rational matrix is evaluated on a (regular) grid 
#' of points on the unit circle ($|z|=1$). The choice \code{which = 'real'} is of course designed 
#' for the case that the z's are real numbers.
#' 
#' \code{NA} values for the (optional) parameters mean that the methods use some suitable default values. 
#' (For the frequency response plots, e.g. the labels for the x- and y-axis are chosen according to the 
#' parameter \code{which}.)
#' 
#' \code{NUll} values for the (optional) parameters mean that the respective graphic element is omitted.
#' E.g. \code{subfigure_main=NULL} skips the titles for the subfigures. 
#' 
#' The titles for the (m,n) subfigures are determined by the parameter \code{subfigure_main}. 
#' One may pass an \eqn{(m,n)} matrix of character strings to the procedure. Alternatively one may also 
#' provide an expression vector with a \code{dim} attribute such that {subfigure_main[i,j]} 
#' returns the expression which is to be used as title for the \eqn{(i,j)}-th subfigure. 
#' E.g. for a $(2,2)$ rational matrix one might use 
#' \preformatted{
#' subfigure_main = expression(Alpha, Beta, Gamma, Delta)
#' dim(subfigure_main) = c(2,2)
#' }
#' 
#' If \code{subfigure_main} is a scalar (character string or expression) then the procedures creates 
#' the respective titles by replacing the "place holders" \code{i_} and \code{j_} with the 
#' respective row and column number. See the examples below.
#' 
#' The parameter \code{xlim} determines the x-axis of the subfigures. 
#' If \code{xlim='subfigure'} then each subfigure gets its 
#' own x-axis (with different limits). For the case \code{xlim='column'} the subfigues in  
#' each column share a common x-axis (with common limits). The case \code{xlim = 'global'} 
#' means that all subfigures have the (x-) limits and the x-axis is only displayed 
#' for the subfigures in the last row. The parameter \code{ylim} is handled analogously. 
#' 
#' If more than one object is plotted (if the optional parameter \code{x_list} is not empty)
#' then a suitable legend may be added with the parameters \code{legend, legend_args}. 
#' Suppose that \eqn{k} objects are plotted, then \code{legend} should be a 
#' character (or expression) vector of length \code{k}. 
#' 
#' The "style" parameters \code{col, type, ..., bg.points} determine the appearance of the 
#' "lines" for each of the \code{k} objects. (If necessary the values are "recycled".)
#' 
#' These plot methods use the internal helper function \code{\link{plot_3D}}. 
#'
#' @param x \code{\link{pseries}} or \code{\link{zvalues}} object.
#' @param x_list (optional) list of additional \code{pseries} or \code{zvalues} objects.
#' @param which (character string) what to plot (only used for \code{\link{zvalues}} 
#'              objects). See details.
#' @param main (character or \code{\link{expression}}) main title of the plot
#' @param ... not used.
#' @inheritParams plot_3D
#'
#' @return The plot methods return (invisibly) a function, \code{subfig} say, which may be used to 
#'         add additional graphic elements to the subfigures. The call \code{opar = subfig(i,j)}
#'         creates a new (sub) plot at the (i,j)-the position with suitable margins and 
#'         axis limits. See the examples below.
#' @export
#'
#' @rdname plot
#' @name plot methods
plot.pseries = function(x, x_list = NULL, 
                        xlim = c('global','column','subfig'), ylim = c('row','subfig','global'),
                        main = 'impulse response', xlab = 'lag (k)', ylab = NULL, 
                        subfigure_main = NA, parse_subfigure_main = FALSE,
                        style = c('gray', 'bw', 'bw2', 'colored'), 
                        col = NA, type = 'l', lty = 'solid', lwd = 1, 
                        pch = 16, cex.points = 1, bg.points = 'black',
                        legend = NULL, legend_args = NA, ...) {
  style = match.arg(style)
  xlim = match.arg(xlim)
  ylim = match.arg(ylim)
  
  ir = unclass(x)
  m0 = dim(ir)[1]
  n0 = dim(ir)[2]
  lag.max  = dim(ir)[3] - 1
  if ((n0*m0*(lag.max+1)) == 0) stop('empty impulse response')
  
  y = list(ir)
  x = list(0:lag.max)
  
  if ((!is.null(x_list)) && (is.list(x_list))) {
    for (i in (1:length(x_list))) {
      if (!is.pseries(x_list[[i]])) stop('argument "x_list" must be a list of "pseries" objects')
      ir = unclass(x_list[[i]])
      m = dim(ir)[1]
      n = dim(ir)[2]
      lag.max  = dim(ir)[3] - 1
      if ( (m != m0) || (n != n0)) stop('impulse response objects are not compatible')
      if ((lag.max+1) == 0) stop('empty impulse response')
      
      y = c(y, list(ir))
      x = c(x, list(0:lag.max))
    }
  }
  k = length(x)
  
  if ((!is.null(subfigure_main)) && (!is.expression(subfigure_main)) && is.na(subfigure_main)) {
    if ((m0*n0) > 1) {
      subfigure_main = '(i_,j_)-th entry'
    } else {
      subfigure_main = NULL
    }
  }
  
  subfigure = plot_3D(x, y, 
                      xlim = xlim, ylim = ylim,  
                      main = main, ylab = ylab, xlab = xlab, 
                      subfigure_main = subfigure_main, parse_subfigure_main = parse_subfigure_main,
                      style = style, 
                      col = col, type = type, lty = lty, lwd = lwd, pch = pch,
                      legend = legend, legend_args = legend_args)  
  return(invisible(subfigure))
}




#' @rdname plot
#' @export
plot.zvalues = function(x, x_list = NULL,  style = c('gray', 'bw', 'colored'), 
                        which = c('modulus','phase','nyquist','real'),
                        subfigure_main = NA, 
                        xlim = NA, ylim = NA, 
                        main = NA, ylab = NA, xlab = NA, 
                        legend = NULL, legend_args = NA, 
                        col = NA, type = 'l', lty = 'solid', lwd = 1, 
                        pch = 16, cex.points = 1, bg.points = 'black', ...) {
  style = match.arg(style)
  which = match.arg(which)
  
  z = attr(x, 'z')
  n.z = length(z)
  
  fr = unclass(x)
  attr(fr, 'z') = NULL
  m0 = dim(fr)[1]
  n0 = dim(fr)[2]
  if ((n0*m0*n.z) == 0) stop('empty frequency response')
  
  x = list(z)
  y = list(fr)
  
  if ((!is.null(x_list)) && (is.list(x_list))) {
    for (i in (1:length(x_list))) {
      if (!is.zvalues(x_list[[i]])) stop('argument "x_list" must be a list of "zvalues" objects')
      
      z = attr(x_list[[i]],'z')
      n.z = length(z)
      
      fr = unclass(x_list[[i]])
      attr(fr, 'z') = NULL
      m = dim(fr)[1]
      n = dim(fr)[2]
      if ( (m != m0) || (n != n0)) stop('zvalues objects are not compatible')
      if (n.z == 0) stop('empty frequency response')
      
      x = c(x, list(z))
      y = c(y, list(fr))
    }
  }
  k = length(x)
  
  if (which == 'modulus') {
    for (i in (1:k)) {
      x[[i]] = -Arg(x[[i]])/(2*pi)
      y[[i]] = abs(y[[i]])
      o = order(x[[i]])
      x[[i]] = x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if (is.na(main)) main = 'frequency response'
    if (is.na(xlab)) xlab = 'frequency (f)'
    if (is.na(ylab)) ylab = 'modulus - frequency response'
    if (is.na(xlim)) {
      xlim = 'column'
    } else {
      xlim = match.arg(xlim, c('global','column','subfig'))
    }
    if (is.na(ylim)) {
      ylim = 'row'
    } else {
      ylim = match.arg(ylim, c('global','row','subfig'))
    }
  }
  if (which == 'phase') {
    for (i in (1:k)) {
      x[[i]] = -Arg(x[[i]])/(2*pi)
      y[[i]] = Arg(y[[i]])/(2*pi)
      o = order(x[[i]])
      x[[i]] = x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if (is.na(main)) main = 'frequency response'
    if (is.na(xlab)) xlab = 'frequency (f)'
    if (is.na(ylab)) ylab = 'argument - frequency response'
    if (is.na(xlim)) {
      xlim = 'column'
    } else {
      xlim = match.arg(xlim, c('global','column','subfig'))
    }
    if (is.na(ylim)) {
      ylim = 'row'
    } else {
      ylim = match.arg(ylim, c('global','row','subfig'))
    }
  }
  if (which == 'nyquist') {
    for (i in (1:k)) {
      d = dim(y[[i]])
      # "close" the curves
      x[[i]] = Re(y[[i]])[,,c(1:d[3],1), drop = FALSE]
      y[[i]] = Im(y[[i]])[,,c(1:d[3],1), drop = FALSE]
    }
    if (is.na(main)) main = 'Nyquist plot'
    if (is.na(xlab)) xlab = 'real part - frequency response'
    if (is.na(ylab)) ylab = 'imaginary part - frequency response'
    if (is.na(xlim)) {
      xlim = 'subfig'
    } else {
      xlim = match.arg(xlim, c('global','column','subfig'))
    }
    if (is.na(ylim)) {
      ylim = 'subfig'
    } else {
      ylim = match.arg(ylim, c('global','row','subfig'))
    }
  }
  if (which == 'real') {
    for (i in (1:k)) {
      x[[i]] = Re(x[[i]])
      y[[i]] = Re(y[[i]])
      o = order(x[[i]])
      x[[i]] = x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if (is.na(main)) main = 'plot - function values'
    if (is.na(xlab)) xlab = 'x'
    if (is.na(ylab)) ylab = 'f(x)'
    if (is.na(xlim)) {
      xlim = 'column'
    } else {
      xlim = match.arg(xlim, c('global','column','subfig'))
    }
    if (is.na(ylim)) {
      ylim = 'row'
    } else {
      ylim = match.arg(ylim, c('global','row','subfig'))
    }
  }
  
  if ((!is.null(subfigure_main)) && (!is.expression(subfigure_main)) && is.na(subfigure_main)) {
    if ((m0*n0) > 1) {
      subfigure_main = '(i_,j_)-th entry'
    } else {
      subfigure_main = NULL
    }
  }
  subfigure = plot_3D(x, y, style = style, subfigure_main = subfigure_main, 
                      xlim = xlim, ylim = ylim,  
                      main = main, ylab = ylab, xlab = xlab, legend = legend, legend_args = legend_args, 
                      col = col, type = type, lty = lty, lwd = lwd, pch = pch)  
  return(invisible(subfigure))
}

