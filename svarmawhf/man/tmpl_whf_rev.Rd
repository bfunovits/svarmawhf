% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aa_template_initvals.R
\name{tmpl_whf_rev}
\alias{tmpl_whf_rev}
\title{Generate an ARMA-WHF template}
\usage{
tmpl_whf_rev(
  dim_out,
  ARorder,
  MAorder,
  kappa,
  k,
  shock_distr = c("laplace", "mixture", "gaussian", "sgt", "skewed_t", "skewed_ged")
)
}
\arguments{
\item{dim_out}{Integer. Output dimension.}

\item{ARorder}{Integer. Degree of AR polynomial matrix}

\item{MAorder}{Integer. Degree of MA polynomial matrix}

\item{kappa, k}{Integers. The partial indices are such that \eqn{\left(\kappa + 1, \ldots, \kappa+1, \kappa, \ldots, \kappa \right)}, 
where the first \code{k} partial indices are equal to \eqn{\kappa+1}.}

\item{shock_distr}{Character string. One out of \code{c("laplace", "mixture", "gaussian", "sgt", "skewed_t")}.}
}
\value{
List with slots \code{ar}, \code{ma_bwd}, \code{ma_fwd}, \code{B} which contain in turn the slots
  \itemize{
    \item \code{h}: vector specifying affine shift of the template
    \item \code{H}: matrix in the affine mapping from deep to linear parameters
    \item \code{n_par}: number of free parameters
    \item \code{degree}: polynomial degree of polynomial matrix
  }
  Furthermore, there are slots
  \itemize{
    \item \code{distr}: A list with slots \code{n_par} and \code{shock_distr}. 
      The first slot is an integer equal to the number of parameters necessary to specify the distribution (zero for Gaussian, three for SGT). 
      The second slot is a string as given by argument \code{shock_distr}.
    \item \code{input_orders}: list with slots \code{dim_out}, \code{ARorder}, \code{MAorder}, \code{kappa}, \code{k}.
    \item \code{n_par}: Total number of free parameters (AR, MA backward, MA forward, B, shock distribution)
  }
}
\description{
Called by main script \code{scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R}.
\cr
For given integer valued parameters \code{dim_out}, \code{ARorder}, \code{MAorder}, \code{kappa}, \code{k}, and given \code{shock_distr}, 
a template connecting the deep/free parameters to the linear parameters is created.
\cr
A template is specified by an affine mapping from the deep/free parameters in the \strong{normalised canonical WHF} to the linear parameters, 
described by a matrix and a shift vector \eqn{(H,h)}.
This mapping is specified separately for the AR, backward MA, forward MA, B matrix.
The template can be filled with deep parameters by using \link{fill_tmpl_whf_rev} and 
the deep parameters can be extracted from a SVARMA-WHF (see return value of \link{fill_tmpl_whf_rev}) object by \link{extract_deep_rev}.
\cr
Note that the variance of each distribution is equal to one. 
Therefore, the matrix \code{B} is unrestricted in the estimation and the permuation/scaling will be obtained afterwards.
}
\examples{
dim_out = 2
ARorder = 3
MAorder = 2
kappa = 1
k = 1
shock_distr = "laplace"
out = tmpl_whf_rev(dim_out = dim_out,
                   ARorder = ARorder,
                   MAorder = MAorder,
                   kappa = kappa,
                   k = k,
                   shock_distr = shock_distr)
out
}
