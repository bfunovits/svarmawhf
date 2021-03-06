% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ba_llfactory.R
\name{get_residuals_once}
\alias{get_residuals_once}
\title{Residuals from deep parameters}
\usage{
get_residuals_once(params_deep, tmpl, data_long)
}
\arguments{
\item{params_deep}{vector of doubles. System and noise parameters which fill the SVARMA-WHF template.}

\item{tmpl}{List. Affine mappings etc specifiying a SVARMA-WHF template. See return value of \link{tmpl_whf_rev}.}

\item{data_long}{matrix. data in long format, i.e. every variable corresponds to a column}
}
\value{
Residuals (long format)
}
\description{
Called in \link{create_results_list}.
}
\details{
Different input arguments than \link{get_residuals} which it calls after input transformation.
Used when calculating the residuals for the Laplace (after call to \link{get_ic}) and SGT version for the first time.
}
