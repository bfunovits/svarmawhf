#' Identification and Estimation of Possibly Non-Invertible SVARMA Models: The Normalised Canonical WHF Parametrisation
#'
#' This is the associated R-package to the paper with the above title.
#' 
#' @section Abstract:
#' 
#' This article deals with parametrisation, identifiability, and quasi maximum likelihood (QML) estimation of possibly non-invertible structural vector autoregressive moving average (SVARMA) models driven by independent and non-Gaussian shocks. 
#' In contrast to previous literature, the novel representation of the MA polynomial matrix using the Wiener-Hopf factorisation (WHF) focuses on the multivariate nature of the model, generates insights into its structure, and uses this structure for devising optimisation algorithms. 
#' In particular, the WHF allows to parameterise the number of MA zeros inside the unit circle, including MA zeros at zero which can be interpreted as informational delays. 
#' This is highly relevant for data-driven evaluation of Dynamic Stochastic General Equilibrium (DSGE) models. 
#' Typically imposed identifying restrictions on the shock transmission matrix as well as on the determinantal root location are made testable. 
#' Furthermore, we provide low level conditions for asymptotic normality of the QML estimator and analytic expressions for the score and the information matrix. 
#' As application, we estimate the Blanchard and Quah model and compare our results to previous literature. 
#' This novel methods are implemented in a well-documented R-package.
#' 
#' @section Dependencies:
#' 
#' This package depends heavily on the \strong{rationalmatrices} and the \strong{RLDM} packages developed by Wolfgang Scherrer and Bernd Funovits.
#' Since both packages may be subject to breaking changes, all necessary functionalities are extracted in order to provide stable code with minimal external dependencies.
#' Importantly, please reference the \strong{rationalmatrices} and the \strong{RLDM} packages should you use their functionalities and not this package.
#' 
#' @section Author(s):
#' 
#' Bernd Funovits and Wolfgang Scherrer
#'
#' Maintainer: <bernd.funovits@@gmail.com>
#'
#' @docType package
#' @name svarmawhf

#' @useDynLib svarmawhf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Pipe operator
#'
#' Re-export pipe operator \code{\%>\%} to turn function composition into a series of imperative statements.
#' For more extensive description, see function \code{`\%>\%`} in package \emph{magrittr}.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs First argument of the function of the right-hand-side of the pipe operator.
#' @param rhs Function whose first argument is given by the left-hand-side argument \code{lhs} of the pipe operator.
#'
#' @examples
#' x = array(stats::rnorm(2*1*3, sd = 0.01), dim = c(2,1,3))
#' # Instead of
#' polm(x)
#' # you can write
#' x %>% polm()
NULL
