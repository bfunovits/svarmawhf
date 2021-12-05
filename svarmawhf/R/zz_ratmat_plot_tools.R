# plot_tools.R ##################################################

#' Plot 3D Arrays
#'
#' This internal helper function plots data arranged in 3-dimensional arrays.
#' The base application is as follows. Let x,y be two 3-dimensional arrays
#' with dimension c(m,n,l). Then \code{plot_3D} splits the device region into a
#' a matrix like array of (m times n) subfigures. In the (i,j)-th subfigure
#' \code{y[i,j,]} is plotted against \code{x[i,j,]}.
#'
#' This function is mainly used by the plot methods, see \code{\link{plot methods}}.
#'
#' For more information about the parameters \code{col, type, lty, ...,bg.points}
#' see \code{\link{plot.default}}.
#'
#' This helper function only makes some very basic checks on the given parameters.
#'
#' @param x list of 3D-arrays (define the 'x'-values)
#' @param y list of 3D-arrays (define the 'y'-values)
#' @param xlim,ylim determine the axis limits of the subfigures.
#'          E.g. \code{xlim = 'column'} means that all subfigures in a column
#'          use the same x-axis limits. The parameter \code{xlim} may also
#'          contain a 2-dimensional vector \code{c(x1,x2)}. In this case
#'          all sub-figures use the given limits for the x-axis. Furthermore
#'          the limits for the y-axis are computed based on the corresponding
#'          data "subset".
#' @param log a character string which contains "x" if the x axis is to be logarithmic,
#'            "y" if the y axis is to be logarithmic and "xy" or "yx" if both axes
#'            are to be logarithmic.
#' @param main (character or \code{\link{expression}}) main title of the plot
#' @param xlab (character string or \code{\link{expression}}) label for the x-axis
#' @param ylab (character or \code{\link{expression}}) label for the y-axis
#' @param subfigure_main scalar or \code{(m x n)} matrix of type "character"
#'             with the titles for the subfigures. If subfigure_main is a scalar
#'             character string then the procedures creates a matrix of respective titles
#'             by replacing the "place holders" '\code{i_}' and '\code{j_}'
#'             with the respective row and column number.
#' @param parse_subfigure_main boolean. If \code{TRUE} then the titles for the subfigures
#'            are parsed to \code{\link{expression}} before plotting. See also
#'            \code{\link[grDevices]{plotmath}} on the usage of expressions for
#'            plot annotations.
#' @param style (character string) determines the appearance of the plot
#'            (background color of the plot regions, color and line style of the
#'            grid lines, axis color, ...) See also \code{\link{style_parameters}}.
#' @param col vector of line colors
#' @param type vector of plot types. The following values are possible: 
#'   "p" for points,
#'   "l" for lines, 
#'   "b" for both points and lines, 
#'   "c" for empty points joined by lines,
#'   "o" for overplotted points and lines, 
#'   "s" and "S" for stair steps and
#'   "h" for histogram-like vertical lines. 
#'   \code{'n'} suppresses plotting.
#' @param lty vector of line types. Line types can either be specified as integers
#'        (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash)
#'        or as one of the character strings "blank", "solid", "dashed", "dotted", "dotdash",
#'        "longdash", or "twodash", where "blank" uses ‘invisible lines’ (i.e., does not draw them).
#' @param lwd vector of line widths.
#' @param pch vector of plotting character or symbols. See \code{\link{points}} for possible
#'            values.
#' @param cex.points vector of scales for the plotting symbols.
#' @param bg.points vector of fill color for the open plot symbols.
#' @param legend (character or \code{\link{expression}} vector). If \code{NULL} then
#'        no legend is produced.
#' @param legend_args (optional) list with parameters for the legend. 
#'   A legend title can be included with \code{legend_args = list(title = my_legend_title)}.
#'   Note that the slots \code{x}, \code{y} are ignored and the legend is always put at the
#'   right hand side of the plot. 
#'   See also \code{\link[graphics]{legend}}.
#'
#' @return A \code{\link{function}} ("closure"), \code{subfig} say, which may be to add additional
#'         graphical elements to the subfigures. The call \code{opar = subfig(i,j)} sets up the
#'         coordinate system margins and figure coordinates such that one may add lines, text, points, ...
#'         to the \code{(i,j)}-th subfigure.
#'
#' @export
plot_3D = function(x, y,
                   xlim = c('subfig','column','global'), ylim = c('subfig','row','global'), log = '',
                   main = NULL, xlab = NULL, ylab = NULL,
                   subfigure_main = '(i_,j_)-th entry', parse_subfigure_main = FALSE,
                   style = c('gray', 'bw', 'bw2', 'colored'),
                   col = NA, type = 'l', lty = 'solid', lwd = 1,
                   pch = 16, cex.points = 1, bg.points = 'black',
                   legend = NULL, legend_args = NA) {
  style = match.arg(style)
  
  # __Call to style_parameters() ####
  style = style_parameters(style)
  
  if (!is.list(x)) x = list(x)
  k_x = length(x)
  if (!is.list(y)) y = list(y)
  k_y = length(y)
  
  if ((k_x == 0) || (k_y == 0)) stop('no x or no y data supplied')
  if (!( (k_x == 1) || (k_y == 1) || (k_x == k_y) )) stop('the lists "x","y" are not compatible')
  k = max(k_x, k_y)
  if ((k_x == 1) && (k > 1)) x = x[rep(1, k)]
  if ((k_y == 1) && (k > 1)) y = y[rep(1, k)]
  
  # only minimal input checks!
  chk = function(x) {
    (is.numeric(x) || is.complex(x)) && (length(x) > 0) && ( (is.vector(x)) || (is.array(x) && (length(dim(x)) == 3)) )
  }
  if (!all(sapply(x, FUN = chk))) {
    stop('argument "x" must be a list of (non-empty, numeric or complex) vectors or 3D arrays')
  }
  if (!all(sapply(y, FUN = chk))) {
    stop('argument "y" must be a list of (non-empty, numeric or complex) vectors or 3D arrays')
  }
  
  # get dimensions (m,n)
  dim2 = function(x) {
    if (is.vector(x)) return(c(1,1))
    return(dim(x)[1:2])
  }
  dim_x = sapply(x, FUN = dim2)
  dim_y = sapply(y, FUN = dim2)
  junk = apply(cbind(dim_x, dim_y), FUN = max, MARGIN = 1)
  m = junk[1]
  n = junk[2]
  
  # convert vectors to 3-D arrays
  conv3D = function(x, m, n) {
    if (is.vector(x)) {
      x = array(x, dim = c(length(x), m, n))
      x = aperm(x, c(2,3,1))
    }
    d = dim(x)[1:2]
    if (!all(d == c(m,n))) stop('array/vector is not compatible')
    return(x)
  }
  x = try(lapply(x, FUN = conv3D, m, n), silent = TRUE)
  if (inherits(x, 'try-error')) stop('x arrays are not compatible')
  y = try(lapply(y, FUN = conv3D, m, n), silent = TRUE)
  if (inherits(y, 'try-error')) stop('y arrays are not compatible')
  
  # check dim[3]
  for (i in (1:k)) {
    l_x = dim(x[[min(i, k_x)]])[3]
    l_y = dim(y[[min(i, k_y)]])[3]
    if ( l_x != l_y ) stop('x, y arrays are not compatible')
  }
  
  # data range
  #range_3D = function(x, MARGIN = integer(0)) {
  X = do.call(dbind, c(list(d=3), x))
  if (grepl('x', log)) X[X<=0] = NA_real_
  Y = do.call(dbind, c(list(d=3), y))
  if (grepl('y', log)) Y[Y<=0] = NA_real_
  if (is.character(xlim)) {
    xlim = match.arg(xlim, c('global','column','subfig'))
    
    # __Call to range_3D() ####
    x_range = switch(xlim,
                     global = range_3D(X),
                     column = range_3D(X, MARGIN = 2),
                     subfig = range_3D(X, MARGIN = 1:2))
  } else {
    xlim = as.vector(xlim)
    if ( (!is.numeric(xlim)) || (length(xlim) != 2) || any(!is.finite(xlim)) ) {
      stop('illegal "xlim" argument!')
    }
    if ( grepl('x', log) && any(xlim <= 0) ) stop('"xlim" has non negative entries and logarithmix x-axis')
    xlim = sort(xlim)
    # print(Y[1,1,])
    Y[(X < xlim[1]) | (X > xlim[2])] = NA_real_
    # print(Y[1,1,])
    x_range = aperm(array(xlim, c(2,m,n)), c(2,3,1))
    xlim = 'global'
  }
  ylim = match.arg(ylim, c('global','row','subfig'))
  y_range = switch(ylim,
                   global = range_3D(Y),
                   row = range_3D(Y, MARGIN = 1),
                   subfig = range_3D(Y, MARGIN = 1:2))
  
  # __Call to default_colmap() ####
  # ()||() and ()&&() conditions must have length 1
  if (is.null(col) || all(is.na(col))) col = default_colmap(k)
  col = as.vector(col)
  col = col[(0:(k-1)) %% length(col) + 1]
  
  type = as.vector(type)
  type = type[(0:(k-1)) %% length(type) + 1]
  lty = as.vector(lty)
  lty = lty[(0:(k-1)) %% length(lty) + 1]
  lwd = as.vector(lwd)
  lwd = lwd[(0:(k-1)) %% length(lwd) + 1]
  pch = as.vector(pch)
  pch = pch[(0:(k-1)) %% length(pch) + 1]
  cex.points = as.vector(cex.points)
  cex.points = cex.points[(0:(k-1)) %% length(cex.points) + 1]
  bg.points = as.vector(bg.points)
  bg.points = bg.points[(0:(k-1)) %% length(bg.points) + 1]
  
  # deal with subfigure_main
  if (!is.null(subfigure_main)) {
    
    stopifnot('Input argument *subfigure_main* must be of type *character*' = is.character(subfigure_main)) 
    
    if (length(subfigure_main) == 1) {
      subfigure_main = matrix(subfigure_main, nrow = m, ncol = n)
      
      # replace 'place holder' i_, j_ by i,j
      for (i in (1:m)) {
        for (j in (1:n)) {
          subfigure_main[i,j] = gsub('i_',paste(i), subfigure_main[i,j])
          subfigure_main[i,j] = gsub('j_',paste(j), subfigure_main[i,j])
        }
      }
    }
    
    stopifnot("Input argument *subfigure_main* must be a matrix (or scalar of type character) with appropriate dimensions." = 
                is.matrix(subfigure_main) && all(dim(subfigure_main) == c(m, n)))
  }
  
  if (!is.null(legend)) {
    if (!is.list(legend_args)) {
      legend_args = list(legend = legend, fill = col, border = col, bty = 'n')
    }
    legend_args$legend = legend
    # ()||() and ()&&() conditions must have length 1
    if ( (!is.null(legend_args$fill)) && all(is.na(legend_args$fill)) ) legend_args$fill = col
    if ( (!is.null(legend_args$border)) && all(is.na(legend_args$border)) ) legend_args$border = col
    if ( (!is.null(legend_args$col)) && all(is.na(legend_args$col)) ) legend_args$col = col
    if ( (!is.null(legend_args$lty)) && all(is.na(legend_args$lty)) ) legend_args$lty = lty
    if ( (!is.null(legend_args$lwd)) && all(is.na(legend_args$lwd)) ) legend_args$lwd = lwd
    if ( (!is.null(legend_args$pch)) && all(is.na(legend_args$pch)) ) legend_args$pch = pch
    if ( (!is.null(legend_args$pt.cex)) && all(is.na(legend_args$pt.cex)) ) legend_args$pt.cex = cex.points
    if ( (!is.null(legend_args$pt.bg)) && all(is.na(legend_args$pt.bg)) ) legend_args$pt.bg = bg.points
  } else {
    legend_args = NULL
  }
  
  tick_labels = array(FALSE, dim = c(m,n,4))
  if (xlim == 'subfig') {
    tick_labels[,,1] = TRUE
  } else {
    tick_labels[m,,1] = TRUE
  }
  if (ylim == 'subfig') {
    tick_labels[,,2] = TRUE
  } else {
    tick_labels[,1,2] = TRUE
  }
  axes = tick_labels
  axes[,,1] = TRUE
  axes[,,2] = TRUE
  titles = array(NA_character_, dim = c(m,n,4))
  if (!is.null(subfigure_main)) titles[,,3] = subfigure_main
  
  # cat('plot_3D:\n', xlim, ylim, '\n')
  # print(tick_labels[,,1])
  # print(tick_labels[,,2])
  
  # compute margins
  margins = list(bottom = 0.2, left = 0.2, top = 0.2, right = 0.2)
  for (i in (1:4)) {
    margins[[i]] = margins[[i]] +
      pmax(1.2 * apply(matrix(tick_labels[,,i], nrow = m, ncol = n), MARGIN = (i-1) %% 2 + 1, FUN = any),
           0.2 * apply(matrix(axes[,,i], nrow = m, ncol = n), MARGIN = (i-1) %% 2 + 1, FUN = any),
           1 * apply(matrix(!is.na(titles[,,i]), nrow = m, ncol = n), MARGIN = (i-1) %% 2 + 1, FUN = any))
  }
  
  # __Call to set_default_par() ####
  # set default graphic parameters
  opar = set_default_par(m,n)
  
  # __Call to start_plot() ####
  # start the plot
  start_plot(xlab = xlab, ylab = ylab, main = main, legend_args = legend_args)
  
  # __Call to get_subfigure_layout() ####
  # compute the corresponding subfigures layout
  subfigure_layout = get_subfigure_layout(margins)
  
  # print(x_range)
  # print(y_range)
  subfigure_fun = function(i = 1, j = 1) {
    # cat(i,j,'\n')
    # print(str(subfigure_layout))
    opar = graphics::par(omi = subfigure_layout$omi, mar = subfigure_layout$mar[i,j,],
                         fig = subfigure_layout$fig[i,j,], new = TRUE)
    # starting a new plot with plot.window() does not work ????
    graphics::plot(x_range[i,j,], y_range[i,j,], type = 'n', log = log,
                   axes = FALSE, frame.plot = FALSE,
                   xlab = NA, ylab = NA, main = NA, sub = NA)
    # cat('subfig_fun:', par('oma'), '\n', mar, '\n', fig, '\n', xlim, '\n', ylim, '\n')
    # graphics::box(which = 'figure', col = 'red')
    return(invisible(opar))
  }
  
  for (i in (1:m)) {
    for (j in (1:n)) {
      subfigure_fun(i,j)
      
      # __Call to plot_axes() ####
      plot_axes(axes = axes[i,j,], tick_labels = tick_labels[i,j,], x_date = 0,
                titles = titles[i,j,], style = style, parse_titles = parse_subfigure_main)
      
      for (l in (1:k)) graphics::lines(x[[l]][i,j,], y[[l]][i,j,],
                                       type = type[l], col = col[l], lty = lty[l],
                                       lwd = lwd[l], pch = pch[l], cex = cex.points[l],
                                       bg = bg.points[l])
    }
  }
  
  # reset graphic parameters
  graphics::par(opar)
  return(invisible(subfigure_fun))
}

# Level 1 functions ####
#' Style Parameters
#'
#' The plotting functions use a list of common style parameters.
#'
#' @param style character indicating the desired "style".
#'
#' @return A list with slots:
#'   \item{border.grid}{background color of the plot region}
#'   \item{col.grid,lty.grid,lwd.grid}{color, line style and line width of the grid lines}
#'   \item{col.axis}{axis color}
#'   \item{col.box,lwd.box}{color and linew idth  of the box surrounding the plot region}
#'   \item{bg.labels,border.labels}{background and border color for axis labels}
#'   \item{col.labels}{text color of the axis labels}
#'
#' @export
#' @keywords internal
style_parameters = function(style = c('bw', 'bw2', 'gray', 'colored')) {
  style = match.arg(style)
  
  if (style == 'bw') {
    style = list(
      bg.grid = NA,
      border.grid = NA,
      col.grid = grDevices::gray(0.5),
      lty.grid = 'dotted',
      lwd.grid = 1,
      col.axis = 'black',
      col.box = 'black',
      lwd.box = 1,
      bg.labels = NA,
      border.labels = NA,
      col.labels = 'black'
    )
    return(style)
  }
  if (style == 'bw2') {
    style = list(
      bg.grid = NA,
      border.grid = NA,
      col.grid = grDevices::gray(0.5),
      lty.grid = 'dotted',
      lwd.grid = 1,
      col.axis = 'black',
      col.box = 'black',
      lwd.box = 1,
      bg.labels = grDevices::gray(0.75),
      border.labels = 'black',
      col.labels = 'black'
    )
    return(style)
  }
  if (style == 'gray') {
    style = list(
      bg.grid = grDevices::gray(0.9),
      border.grid = NA,
      col.grid = 'white',
      lty.grid = 'solid',
      lwd.grid = 1,
      col.axis = grDevices::gray(0.5),
      col.box = NA,
      lwd.box = 1,
      bg.labels = grDevices::gray(0.5),
      border.labels = NA,
      col.labels = 'white'
    )
    return(style)
  }
  if (style == 'colored') {
    style = list(
      bg.grid = grDevices::gray(0.9),
      border.grid = NA,
      col.grid = 'white',
      lty.grid = 'solid',
      lwd.grid = 1,
      col.axis = grDevices::gray(0.5),
      col.box = NA,
      lwd.box = 1,
      bg.labels = 'orange',
      border.labels = NA,
      col.labels = 'black'
    )
    return(style)
  }
  stop('this should not happen')
}

#' Range of Values in a 3-D Array
#'
#' \code{range_3D} computes the range of values of a three-dimensional array.
#'
#' If the input, \code{x} say, is an \code{(m,n,k)} dimensional array then \code{range_3D}
#' returns an \code{(m,n,2)} dimensional array, \code{m} say:
#'
#' \itemize{
#' \item \code{MARGIN = c(1,2)}: \code{m[i,j,]} contains the range of values in \code{x[i,j,]}.
#' \item \code{MARGIN = 1}: \code{m[i,j,]} contains the range of values in \code{x[i,,]}.
#' \item \code{MARGIN = 2}: \code{m[i,j,]} contains the range of values in \code{x[,j,]}.
#' \item \code{MARGIN = integer(0)}: \code{m[i,j,]} contains the range of values in \code{x[,,]}.
#' }
#'
#' @param x \code{(m,n,k)}-dimensional numeric array
#' @param MARGIN integer(0), 1, 2 or (1:2). determines the "margin" over which to compute the
#'               range(s).
#'
#' @return \code{(m,n,2)}-dimensional array.
#' @export
#' @keywords internal
range_3D = function(x, MARGIN = integer(0)) {
  m = dim(x)[1]
  n = dim(x)[2]
  if (length(MARGIN) == 0) {
    r = suppressWarnings(range(x, na.rm = TRUE)) # 2
    r = array(r, dim = c(2,m,n))
    r = aperm(r, c(2,3,1))
    return(r)
  }
  if ((length(MARGIN) == 1) && (MARGIN == 1)) {
    r = apply(x, MARGIN = 1, FUN = function(m) {suppressWarnings(range(m, na.rm = TRUE))} ) # 2 x m
    r = array(r, dim = c(2,m,n))
    r = aperm(r, c(2,3,1))
    return(r)
  }
  if ((length(MARGIN) == 1) && (MARGIN == 2)) {
    r = apply(x, MARGIN = 2, FUN = function(m) {suppressWarnings(range(m, na.rm = TRUE))} ) # 2 x n
    r = array(r, dim = c(2,n,m))
    r = aperm(r, c(3,2,1))
    return(r)
  }
  if ((length(MARGIN) == 2) && all(MARGIN == c(1,2))) {
    r = apply(x, MARGIN = c(1,2), FUN = function(m) {suppressWarnings(range(m, na.rm = TRUE))} ) # 2 x m x n
    r = aperm(r, c(2,3,1))
    return(r)
  }
  stop('illegal MARGIN parameter')
}


#' Default Color Palette
#'
#' Default color palette.
#'
#' @param n integer
#'
#' @return (character) vector of length \code{n} with color codes (\code{"#RRGGBB"}).
#' @export
#' @keywords internal
#'
#' @examples
#' default_colmap(3)
default_colmap = function(n) {
  return( scales::hue_pal()(n) )
}



#' Set Default Graphical Parameters
#'
#' Set default graphical parameters.
#'
#' @param m,n integers. The procedure assumes that the device is
#'            split up into an m-by-n array of (sub-) figures and scales
#'            \code{cex.axis}, \code{cex.lab} and \code{cex.main}
#'            correspondingy.
#'
#' @return invisible list with the original graphical parameters.
#' @export
#' @keywords internal
#'
#' @examples
#' opar = set_default_par(2,3)
#' print(opar)
set_default_par = function(m=1,n=1){
  sc = sqrt((m^2+n^2)/2)
  sc = max(0.5, 0.5 + 0.5/sqrt(sc))
  opar = graphics::par(oma = c(0, 0, 0, 0),
                       fig = c(0,1,0,1), mar = c(0,0,0,0),
                       tcl = -0.1*sc, mgp = c(1, 0.25, 0)*sc, mex = 1,
                       cex = 1, mex = 1,
                       cex.axis = sc, cex.lab = sc, cex.main = sc,
                       xaxs = 'r', yaxs = 'r', col.axis = 'black')
  return(invisible(opar))
}


#' Start a new Plot
#'
#' This tools starts a new plot and writes the optional "labels"
#' \code{xlab}, \code{ylab} and \code{main} into the (outer) margins
#' of the plot. In addition an (optional) legend is put into the
#' right (outer) margin. The outer margins of the plot are determined
#' based on whether or not text or legend is put into the respective margin.
#'
#' @param xlab,ylab,main (optional) character or \code{\link{expression}}.
#'        If \code{NULL} then the respective outer margin is set to zero.
#' @param legend_args  (optional) list with \code{\link[graphics]{legend}}
#'        arguments. Note that the legend is always put at the "right" side of the
#'        plot. If \code{NULL} then the right outer margin is set to zero.
#'
#' @return (invisible) 4-dimensional vector with the outer margins
#'         \code{par("mar")}.
#' @export
#' @keywords internal
#'
#' @seealso \code{\link[graphics]{legend}} and \code{\link[graphics]{par}}.
#'
#' @examples
#' \dontrun{
#' set_default_par()
#' start_plot(xlab = 'x', ylab = 'y', main = 'main',
#'            legend_args = list(legend = c('eins','zwei'), fill = c('red','blue')))
#'            graphics::box(which = 'inner', col = 'red')
#' }
start_plot = function(xlab = NULL, ylab = NULL,
                      main = NULL, legend_args= NULL) {
  
  # set outer margins (place for xlab, ylab, main)
  opar = graphics::par(xaxs = 'i', yaxs = 'i')
  graphics::par(oma = c((!is.null(xlab))*1.5, (!is.null(ylab))*1.5,
                        (!is.null(main))*1.5*1.2, 0),
                fig = c(0,1,0,1), mar = c(0,0,0,0))
  omd = graphics::par('omd')
  
  # create a new figure which fills the region inside the outer margins
  graphics::plot.new()
  graphics::plot.window(c(0,1), c(0,1))
  
  if (!is.null(legend_args)) {
    # plot legend
    legend_args$x = 'right'
    legend_args$y = NULL
    out = do.call(graphics::legend, legend_args)
    # set right outer margin ( place for the legend )
    omd[2] = out$rect$left
    graphics::par(omd = omd)
  }
  
  # reset axis-styles
  graphics::par(opar)
  
  # plot main, xlab, ylab
  if (!is.null(xlab)) graphics::mtext(xlab, cex = 1,
                                      side = 1, line = 0, outer = TRUE)
  if (!is.null(ylab)) graphics::mtext(ylab, cex = 1,
                                      side = 2, line = 0, outer = TRUE)
  if (!is.null(main)) graphics::mtext(main, cex = 1.2,
                                      side = 3, line = 0.15, outer = TRUE)
  
  return(invisible(graphics::par('mar')))
}





#' Subfigure Layout
#'
#' The device is split up into an (m-by-n) array of (sub-) figures. This tool
#' computes the corresponding figure region and margins for each of the
#' (sub-) figures.
#'
#' @param margins A named list with slots \code{bottom}, \code{left}, \code{top}
#'                and \code{right} which determine the margins of the sub-figures.
#'                E.g. \code{bottom} is an \code{m}-dimensional vector where
#'                \code{bottom[i]} gives the bottom margin (in "line units")
#'                of the sub-figures in the \code{i}-th row.
#'
#' @return A list with slots
#' \item{omi}{4-dimensional vector with the "outer margin" in inches.}
#' \item{mar}{(m,n,4)-dimensional array, where \code{mar[i,j,]} contains the
#'            margins (in "line" units) of the (i,j)-th sub-figure.}
#' \item{fig}{(m,n,4)-dimensional array, where \code{fig[i,j,]} contains the
#'            coordinates (in "NDC" units) of the figure region of the
#'            (i,j)-th sub-figure.}
#'
#' @seealso \code{\link[graphics]{par}} for more details.
#'
#' @export
#' @keywords internal
get_subfigure_layout = function(margins) {
  bottom = rev(margins$bottom)
  top = rev(margins$top)
  left = margins$left
  right = margins$right
  
  dev_size = grDevices::dev.size('in') # width, height of device in inches
  
  omi = graphics::par('omi')  # outer margin in 'inches' and lines
  # cat('dev_size:', dev_size, '\n')
  l2in = line2inch()
  
  # convert margins from lines to inches
  bottom_in = bottom * l2in
  left_in = left * l2in
  top_in = top * l2in
  right_in = right * l2in
  m = length(bottom) # number of rows
  n = length(left)   # number of columns
  plot_width = ( dev_size[1] -
                   (omi[2] + sum(left_in) + sum(right_in) + omi[4]) ) / n
  plot_height = ( dev_size[2] -
                    (omi[1] + sum(bottom_in) + sum(top_in) + omi[3]) ) / m
  # cat('plot_size:', c(plot_width, plot_height), '\n')
  # x coordinate of the lower left corner of the sub figures
  x = cumsum(c(0, left_in + plot_width + right_in))
  # y coordinate of the lower left corner of the sub figures
  y = cumsum(c(0, bottom_in + plot_height + top_in))
  # cat('x:', x, '\n')
  # cat('diff(x):', diff(x),'\n')
  # cat('y:', y, '\n')
  # cat('diff(y):', diff(y),'\n')
  
  # convert to NDC coordinates
  x = x / (dev_size[1] - omi[2] - omi[4])
  x[x > 1] = 1
  y = y / (dev_size[2] - omi[1] - omi[3])
  y[y > 1] = 1
  
  mar = array(0, dim = c(m,n,4))
  fig = array(0, dim = c(m,n,4))
  
  for (i in (1:m)) {
    for (j in (1:n)) {
      mar[i,j, ] = c(bottom[i], left[j], top[i], right[j])
      fig[i,j, ] = c(x[j], x[j+1], y[i], y[i+1])
    }
  }
  mar = mar[m:1,,,drop = FALSE]
  fig = fig[m:1,,,drop=FALSE]
  return(list(omi = omi, mar = mar, fig = fig))
}



#' Plot Grid Lines, Axes and Axes Labels
#'
#' This tools plots the grid lines,axes and axes titles/labels. The x-axis may be
#' a date/time axis.
#'
#' @param axes 4-dimensional boolean vector, indicating whether an axis should be put
#'             at the respective side of the plot.
#' @param tick_labels 4-dimensional boolean vector, indicating whether the respective
#'             axis has tick-labels.
#' @param x_date scalar number, used to convert the (numeric) x-axis limits \code{xlim} to
#'             a date/time object. If \code{x_date} is of class \code{POSIXct}
#'             then the limits are converted as \code{as.POSIXct(xlim, origin = x_date)}.
#'             If \code{x_date} is of class \code{Date}
#'             then the limits are converted as \code{as.Date(xlim, origin = x_date)}.
#'             Otherwise no coercion is performed. Finally
#'             \code{graphics::Axis(x = xlim, ...)} is called to do the actual plotting
#'             of the x-axis.
#' @param titles a 4-dimensional character vector. If the respective entry is \code{NA} then
#'             no label/title is put at this side of the plot.
#' @param style list of "style" parameters as returned by \code{\link{style_parameters}}.
#' @param parse_titles if TRUE the titles are coerced to \code{\link{expression}} before
#'              plotting.
#'
#' @return (invisible) \code{NULL}.
#' @export
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' date = seq(as.Date("2019/1/1"), by = "month", length.out = 24)
#' titles = c(NA, 'y','main', NA)
#' plot(date, 1:24, axes = FALSE, xlab = NA, ylab = NA)
#' plot_axes(axes = c(TRUE, FALSE, FALSE, TRUE), titles = titles,
#'           style = style_parameters('gray'), x_date = date[1] - as.numeric(date[1]))
#' points(date, 1:24)
#' }
plot_axes = function(axes = c(TRUE, TRUE, FALSE, FALSE), tick_labels = axes, x_date = 0,
                     titles = rep(NA_character_, 4), style = style_parameters('bw'),
                     parse_titles = FALSE) {
  junk = graphics::par('cex.axis', 'cex.lab', 'usr', 'xlog', 'ylog')
  cex.axis = junk$cex.axis
  cex.lab = junk$cex.lab
  usr = junk$usr
  xlim = usr[1:2]
  if (junk$xlog) xlim = 10^xlim
  ylim = usr[3:4]
  if (junk$ylog) ylim = 10^ylim
  
  # print(style$bg.grid)
  # cat(xlim[1], ylim[1], xlim[2], ylim[2], style$bg.grid, style$border.grid, '\n')
  # graphics::rect(xlim[1], ylim[1], xlim[2], ylim[2],
  #                col = 'orange')
  # lines(xlim, ylim, xpd = NA)
  graphics::rect(xlim[1], ylim[1], xlim[2], ylim[2],
                 col = style$bg.grid, border = style$border.grid)
  
  at_x = NULL
  at_y = NULL
  for (i in (1:4)) {
    if (axes[i]) {
      if (i == 1) {
        a = 0.2
        delta = 0.25*cex.axis - (1-a)*(1-cex.axis)
      } else {
        delta = 0.25*cex.axis
      }
      if ( (i %% 2) == 1 ) {
        if (inherits(x_date, 'POSIXct')) {
          if (junk$xlog) message('logarithmic scale for "POSIXct" class does not work (well)')
          xlim = as.POSIXct(xlim, origin = x_date)
        }
        if (inherits(x_date, 'Date')) {
          if (junk$xlog) message('logarithmic scale for "Date" class does not work (well)')
          xlim = as.Date(xlim, origin = x_date)
        }
        at_x = graphics::Axis(x = xlim,  at = NULL,
                              lwd = 0, lwd.ticks = 1, col.axis = style$col.axis, mgp = c(0, delta, 0),
                              side = i, labels = tick_labels[i])
      } else {
        at_y = graphics::Axis(x = ylim,  at = NULL,
                              lwd = 0, lwd.ticks = 1, col.axis = style$col.axis, mgp = c(0, delta, 0),
                              side = i, labels = tick_labels[i])
      }
    }
  }
  
  if (!is.null(at_x)) graphics::abline(v = at_x, col = style$col.grid, lty = style$lty.grid, lwd = style$lwd.grid)
  if (!is.null(at_y)) graphics::abline(h = at_y, col = style$col.grid, lty = style$lty.grid, lwd = style$lwd.grid)
  
  if (!is.na(style$col.box)) graphics::box(which = 'plot', col = style$col.box, lwd = style$lwd.box)
  
  
  for (i in (1:4)) {
    if (!is.na(titles[i])) {
      btext(titles[i], side = i, cex = cex.lab, size = cex.lab,
            col = style$col.labels, bg = style$bg.labels, border = style$border.labels,
            parse_text = parse_titles)
    }
  }
  
  return(invisible(NULL))
}


# Level 2 functions ####
#' Write "Boxed" Text into the Margins of a Plot.
#'
#' This tool writes text into a "colored" box in one of the margins
#' of the plot region.
#'
#' @param text character or \code{\link{expression}}.
#' @param side on which side of the plot (1=bottom, 2=left, 3=top, 4=right).
#' @param cex character expansion factor.
#' @param col color to use (for the text).
#' @param size determines the size of the box (in "line" units)
#' @param bg,border background and border color of the box.
#' @param parse_text if yes the procedure tries to coerce the
#'                   text into an \code{\link{expression}}.
#'
#' @return (invisible) vector with the "box" corner and center coordinates.
#' @export
#' @keywords internal
#'
#' @seealso \code{\link[graphics]{mtext}} and \code{\link[grDevices]{plotmath}}
#' for details on mathematical annotation.
#'
#' @examples
#' \dontrun{
#' plot(1:10, axes = FALSE)
#' graphics::box()
#' btext('hallo', side = 1)
#' btext('Sigma[11]', side = 2, bg = 'lightblue')
#' btext('Sigma[11]', side = 3, bg = NA, border = 'black', size = 1.2, parse_text = TRUE)
#' btext(expression(Sigma[11]), side = 4, bg = 'orange', border = 'black', cex = 1.5)
#' }
btext = function(text, side = 3, cex = 1, col = 'black',
                 size = cex, bg = 'lightgray', border = 'lightgray',
                 parse_text = FALSE) {
  l2in = line2inch()
  junk = graphics::par('usr','mai','fin','xlog','ylog')
  usr = junk$usr
  
  plt = junk$fin - c(junk$mai[2]+junk$mai[4],
                     junk$mai[1]+junk$mai[3]) # size of plot region in inches
  in2x = (usr[2]-usr[1])/plt[1] # lines to usr (x) coordinates
  in2y = (usr[4]-usr[3])/plt[2] # lines to usr (y) coordinates
  
  
  # cat(junk$fin, plt, usr, in2x, in2y, '\n')
  
  if ((side %% 2) == 1) {
    srt = 0
    size = size * l2in * in2y
    if (side == 1) {
      box = c(usr[1], usr[2], usr[3] - size, usr[3])
    } else {
      box = c(usr[1], usr[2], usr[4], usr[4] + size)
    }
  } else {
    size = size * l2in * in2x
    if (side == 2) {
      srt = 90
      box = c(usr[1] - size, usr[1], usr[3], usr[4])
    } else {
      srt = -90
      box = c(usr[2], usr[2] + size, usr[3], usr[4])
    }
  }
  center = c(box[1] + box[2], box[3] + box[4] ) /2
  if (junk$xlog) {
    box[1:2] = 10^box[1:2]
    center[1] = 10^center[1]
  }
  if (junk$ylog) {
    box[3:4] = 10^box[3:4]
    center[2] = 10^center[2]
  }
  
  if (is.expression(text)) parse_text = FALSE
  if (parse_text) {
    junk = try(parse(text = text), silent = TRUE)
    if (!inherits(junk, 'try-error')) text = junk
  }
  
  graphics::rect(box[1], box[3], box[2], box[4], col = bg, border = border, xpd = TRUE)
  text(center[1], center[2], text, cex = cex, srt = srt,
       adj = c(0.5, 0.5), col = col, xpd = TRUE)
  
  return(invisible(c(box, center)))
}




#' Convert 'line units' to inches
#'
#' @return scalar
#' @export
#' @keywords internal
#'
#' @examples
#' line2inch()
#' par(cex = 2, mex = 2)
#' line2inch()
line2inch = function() {
  opar = graphics::par(mar = c(1,1,1,1))
  x = mean(graphics::par('mai'))
  graphics::par(opar)
  return(x)
}

# Not sure where this is needed ####

#' Zoom and Scroll
#'
#' The utility \code{zoom_plot} creates and runs a "shiny app" which displays two plots:
#' Drawing (and dragging) a "brush" in the "control" plot at the bottom zooms into
#' zooms into the selected x-range of the plot at the top.
#'
#' This utility is based on the \pkg{shiny} package.
#'
#' @param p,p0  two \code{\link{function}} ("closure") objects which produce an R plot
#'        when calling \code{p()}, \code{p(xlim)} and \code{p0()} respectively.
#'        Both plots should use the same range of "x-values".
#' @param title (character) optimal title for the "shiny app window"
#' @param ... not used.
#'
#' @return This function normally does not return; interrupt R to stop the application
#' (usually by pressing Ctrl+C or Esc).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' make_plot_fun = function(x, y, type, col) {
#'    fun = function(xlim = NULL) {
#'    plot(x, y, xlim= xlim, type = type, col = col)
#'    }
#'    return(fun)
#' }
#'
#' p = make_plot_fun(1:10, rnorm(10), type = 'p', col = 'red')
#' p0 = make_plot_fun(1:10, sin(1:10), type = 'l', col = 'black')
#' zoom_plot(p, p0, title = 'test "zoom_plot"')
#' }
zoom_plot <- function(p, p0 = p, title = NULL, ...) {
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    message('"zoom_plot" needs the "shiny" package')
    return(invisible(NULL))
  }
  
  # create shiny - App
  app <- list(
    ui = shiny::fluidPage(
      shiny::titlePanel(title),
      shiny::fluidRow(style = "background-color: white;",
                      shiny::column(width = 12,
                                    shiny::plotOutput('main',height="500px") # main-plot
                      )
      ),
      shiny::hr(),
      shiny::fluidRow(style = "background-color: white;",
                      # control-plot with a "brush"
                      shiny::column(width = 12,
                                    shiny::plotOutput('control', height="150px",
                                                      brush = shiny::brushOpts(id = "zoom", direction = "x"))
                      )
      )
    ),
    server = function(input, output) {
      output$main <- shiny::renderPlot({
        if (is.null(input$zoom)) {
          p()
        } else {
          # new range of x-values
          xlim = c(input$zoom$xmin, input$zoom$xmax)
          p(xlim)
        }
      })
      output$control <- shiny::renderPlot({
        p0()
      })
      # output$plot_brushinfo <- renderPrint({
      #   cat("Brush (debounced):\n")
      #   str(input$zoom)
      # })
    }
  )
  # run the shiny-app
  shiny::runApp(app, ...)
}



