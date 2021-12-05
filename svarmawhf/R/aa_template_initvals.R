#' Generate an ARMA-WHF template
#' 
#' Called by main script \code{scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R}.
#' \cr
#' For given integer valued parameters \code{dim_out}, \code{ARorder}, \code{MAorder}, \code{kappa}, \code{k}, and given \code{shock_distr}, 
#' a template connecting the deep/free parameters to the linear parameters is created.
#' \cr
#' A template is specified by an affine mapping from the deep/free parameters in the \strong{normalised canonical WHF} to the linear parameters, 
#' described by a matrix and a shift vector \eqn{(H,h)}.
#' This mapping is specified separately for the AR, backward MA, forward MA, B matrix.
#' The template can be filled with deep parameters by using \link{fill_tmpl_whf_rev} and 
#' the deep parameters can be extracted from a SVARMA-WHF (see return value of \link{fill_tmpl_whf_rev}) object by \link{extract_deep_rev}.
#' \cr
#' Note that the variance of each distribution is equal to one. 
#' Therefore, the matrix \code{B} is unrestricted in the estimation and the permuation/scaling will be obtained afterwards.
#'
#' @param dim_out Integer. Output dimension.
#' @param ARorder Integer. Degree of AR polynomial matrix
#' @param MAorder Integer. Degree of MA polynomial matrix
#' @param kappa,k Integers. The partial indices are such that \eqn{\left(\kappa + 1, \ldots, \kappa+1, \kappa, \ldots, \kappa \right)}, 
#'   where the first \code{k} partial indices are equal to \eqn{\kappa+1}.
#' @param shock_distr Character string. One out of \code{c("laplace", "mixture", "gaussian", "sgt", "skewed_t")}.
#'
#' @return List with slots \code{ar}, \code{ma_bwd}, \code{ma_fwd}, \code{B} which contain in turn the slots
#'   \itemize{
#'     \item \code{h}: vector specifying affine shift of the template
#'     \item \code{H}: matrix in the affine mapping from deep to linear parameters
#'     \item \code{n_par}: number of free parameters
#'     \item \code{degree}: polynomial degree of polynomial matrix
#'   }
#'   Furthermore, there are slots
#'   \itemize{
#'     \item \code{distr}: A list with slots \code{n_par} and \code{shock_distr}. 
#'       The first slot is an integer equal to the number of parameters necessary to specify the distribution (zero for Gaussian, three for SGT). 
#'       The second slot is a string as given by argument \code{shock_distr}.
#'     \item \code{input_orders}: list with slots \code{dim_out}, \code{ARorder}, \code{MAorder}, \code{kappa}, \code{k}.
#'     \item \code{n_par}: Total number of free parameters (AR, MA backward, MA forward, B, shock distribution)
#'   }
#'   
#' @export
#' 
#' @examples 
#' dim_out = 2
#' ARorder = 3
#' MAorder = 2
#' kappa = 1
#' k = 1
#' shock_distr = "laplace"
#' out = tmpl_whf_rev(dim_out = dim_out,
#'                    ARorder = ARorder,
#'                    MAorder = MAorder,
#'                    kappa = kappa,
#'                    k = k,
#'                    shock_distr = shock_distr)
#' out
tmpl_whf_rev = function(dim_out, 
                        ARorder, MAorder, kappa, k, 
                        shock_distr = c("laplace", "mixture", "gaussian", "sgt", "skewed_t", "skewed_ged")){
  
  # Polynomial degrees
  stopifnot("Inputs *MAorder*, *kappa*, *k* must be such that *MAorder* is larger than or equal to *kappa* (if *k* is zero) and larger than *kappa* if *k* is larger than zero." = MAorder >= kappa + as.integer(k > 0))
  MAorder_bwd = MAorder - kappa
  MAorder_fwd = kappa + as.integer(k > 0)
  
  # WHF system parameter template ####
  
  # AR
  array_ar = array(NA_real_, 
                   dim = c(dim_out, dim_out, 1 + ARorder))
  array_ar[,,1] = diag(dim_out)
  
  # MA bwd
  array_ma_bwd = array(NA_real_, 
                       dim = c(dim_out, dim_out, 1 + MAorder_bwd))
  array_ma_bwd[,,1] = diag(dim_out)
  
  if (MAorder_bwd > 0 && k > 0){
    array_ma_bwd[,(1:k), 1 + MAorder_bwd] = 0
    array_ma_bwd[(k+1):dim_out, (1:k), 1] = NA_real_
    array_ma_bwd[1:k, (k+1):dim_out, 2] = 0
  }
  
  # MA fwd
  array_ma_fwd = array(NA_real_, 
                       dim = c(dim_out, dim_out, 1 + MAorder_fwd))
  array_ma_fwd[,,1] = diag(dim_out)
  
  if (MAorder_fwd > 0 && k > 0){
    array_ma_fwd[(k+1):dim_out, , 1+MAorder_fwd] = 0
  }
  
  # Affine transformations ####
  
  # AR
  h_ar = as.vector(array_ar)
  ix = which(!is.finite(h_ar))
  h_ar[ix] = 0
  n_sys_ar = length(ix)
  
  H_ar = matrix(0, nrow = length(h_ar), ncol = n_sys_ar)
  H_ar[ix,] = diag(n_sys_ar)
  
  # MA backward
  h_ma_bwd = as.vector(array_ma_bwd)
  ix = which(!is.finite(h_ma_bwd))
  h_ma_bwd[ix] = 0
  n_sys_ma_bwd = length(ix)
  
  H_ma_bwd = matrix(0, nrow = length(h_ma_bwd), ncol = n_sys_ma_bwd)
  H_ma_bwd[ix,] = diag(n_sys_ma_bwd)
  
  # MA forward
  h_ma_fwd = as.vector(array_ma_fwd)
  ix = which(!is.finite(h_ma_fwd))
  h_ma_fwd[ix] = 0
  n_sys_ma_fwd = length(ix)
  
  H_ma_fwd = matrix(0, nrow = length(h_ma_fwd), ncol = n_sys_ma_fwd)
  H_ma_fwd[ix,] = diag(n_sys_ma_fwd)
  
  # B matrix ####
  h_B = double(dim_out^2)
  H_B = diag(dim_out^2)
  n_par_B = dim_out^2
  
  # Distribution ####
  n_par_distr = 
    if (shock_distr %in% c("gaussian", "laplace")){
      0
    } else if (shock_distr %in% c("mixture", "sgt")){
      dim_out * 3
    } else if (shock_distr %in% c("skewed_t", "skewed_ged")){
      dim_out * 2
    } else {
      stop("One of *laplace*, *mixture*, *gaussian*, *sgt*, *skewed_t*, *skewed_ged* needs to be specified in argument *shock_distr*.")
    }
  
  return(list(ar = list(h = h_ar,
                        H = H_ar,
                        n_par = n_sys_ar,
                        degree = ARorder),
              ma_bwd = list(h = h_ma_bwd,
                            H = H_ma_bwd,
                            n_par = n_sys_ma_bwd,
                            degree = MAorder_bwd),
              ma_fwd = list(h = h_ma_fwd,
                            H = H_ma_fwd,
                            n_par = n_sys_ma_fwd,
                            degree = MAorder_fwd),
              B = list(h = h_B,
                       H = H_B,
                       n_par = n_par_B),
              distr = list(n_par = n_par_distr,
                           shock_distr = shock_distr),
              input_orders = list(dim_out = dim_out,
                                  ARorder = ARorder, 
                                  MAorder = MAorder, 
                                  kappa = kappa, 
                                  k = k),
              n_par = n_sys_ar + n_sys_ma_bwd + n_sys_ma_fwd + n_par_B + n_par_distr))
}

#' Obtain Initial deep parameters
#' 
#' Called by main script \code{scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R}.
#' Calls \link{est_ar}, \link{extract_deep_rev}, \link{get_residuals}.
#' \cr
#' Obtains initial deep parameters for given data \code{data_long} and a given SVARMAWHF template \code{tmpl} obtained from \link{tmpl_whf_rev}.
#' The AR parameters are obtained from the Yule-Walker equations, see \link{est_ar}. 
#' In order to ensure that the initial value for the backward and forward MA polynomial matrix be stable (have its determinantal roots on the correct side of the unit circle),
#' we initial them with the (possibly randomized) Minnesota prior \link[LaplacesDemon]{MinnesotaPrior}.
#' \cr
#' Subsequently, we obtain the residuals for the template filled with the initial deep parameter values via \link{get_residuals}.
#'
#' @section Note:
#' Case \code{q = p} throws an error because otherwise problems with \link[stats]{lsfit} in \link{extract_deep_rev} would occur. 
#' Handling in the code would be nicer than in the Rmarkdown file...
#' 
#' @param data_long Matrix. Contains observation in each row, number of variables corresponds to number of rows. 
#'   Note that the data need to be \strong{demeaned}. This is not checked!
#' @param tmpl SVARMAWHF template. As obtained from \link{tmpl_whf_rev}.
#'
#' @return Vector of deep system and noise parameters. 
#'   They are used to fill a given template.
#'   
#' @export
get_init_armamod_whf_random = function(data_long, tmpl){
  
  ARorder = tmpl$ar$degree
  MAorder_bwd = tmpl$ma_bwd$degree
  MAorder_fwd = tmpl$ma_fwd$degree
  
  dim_out = tmpl$input_orders$dim_out
  n_obs = dim(data_long)[1] 
  kappa = tmpl$input_orders$kappa
  k = tmpl$input_orders$k
  
  lambda_this = 0.8
  is_unstable = TRUE
  
  while(is_unstable && lambda_this > 0){
    
    lambda_this = lambda_this - 0.05
    stopifnot("get_init_armamod_whf_random(): Initial values of MA (fwd or bwd) poly are not stable." = lambda_this > 0)
    
    # Stable MA part "past" ####
    polm_ma_bwd = LaplacesDemon::MinnesotaPrior(J = dim_out, 
                                                lags = 1:MAorder_bwd, 
                                                lambda = lambda_this, theta = 0.2, sigma = 1)
    polm_ma_bwd = array(c(c(diag(dim_out)),
                          -c(polm_ma_bwd)),
                        dim = c(dim_out, dim_out, MAorder_bwd+1))
    
    if (k > 0){
      polm_ma_bwd[(k+1):dim_out, 1:k, 1] = c(polm_ma_bwd[1:k,(k+1):dim_out,2])
      polm_ma_bwd[1:k, (k+1):dim_out, 2] = 0
      polm_ma_bwd[, 1:k, 1+MAorder_bwd] = 0
    }
    
    if (MAorder_bwd > 0){
      is_unstable_bwd = any(abs(eigen(companion_matrix(polm_ma_bwd), only.values = TRUE)$values) > 1)
    } else {
      is_unstable_bwd = FALSE
    }
    
    # "Unstable" MA part "future" ####
    polm_ma_fwd = LaplacesDemon::MinnesotaPrior(J = dim_out, 
                                                lags = 1:MAorder_fwd, 
                                                lambda = lambda_this, theta = 0.2, sigma = 1)
    polm_ma_fwd = array(c(c(diag(dim_out)),
                          -polm_ma_fwd),
                        dim = c(dim_out, dim_out, MAorder_fwd + 1))
    if (k > 0){
      polm_ma_fwd[(k+1):dim_out, , 1+MAorder_fwd] = 0
    }
    
    if (MAorder_fwd > 0){
      is_unstable_fwd = any(abs(eigen(companion_matrix(polm_ma_fwd), only.values = TRUE)$values) > 1)
    } else {
      is_unstable_fwd = FALSE
    }
    
    is_unstable = is_unstable_bwd || is_unstable_fwd
  }
  
  
  # AR and "noise" parameters####
  ar_yw_obj = est_ar(data_long,
                     p.max = ARorder,
                     ic = "max",
                     method = "yule-walker", 
                     mean_estimate = "zero")
  # ...$sys holds the lmfd-object in wide-matrix form. 
  # The MA coefficients are implicitly discarded. Not really a safe way...
  polm_ar = ar_yw_obj$model$sys %>% array(dim = c(dim_out, dim_out, ARorder + 1))
  
  if (ARorder > 0){
    stopifnot("get_init_armamod_whf_random(): AR polynomial matrix not stable." = all(abs(eigen(companion_matrix(polm_ar), only.values = TRUE)$values) < 1))
  }
  
  param_sys = extract_deep_rev(armamod_whf_obj = list(polm_ar = polm_ar,
                                                      polm_ma_bwd = polm_ma_bwd,
                                                      polm_ma_fwd = polm_ma_fwd),
                               tmpl = tmpl)
  
  # Noise B: Input preparation for RcppArmadillo function get_residuals() ####
  polm_ar4cpp = matrix(polm_ar, 
                       nrow = dim_out, ncol = dim_out*(1 + tmpl$ar$degree))
  polm_ma_bwd4cpp = matrix(polm_ma_bwd[,, rev(iseq(1, 1+tmpl$ma_bwd$degree))], 
                           nrow = dim_out, ncol = dim_out*(1+tmpl$ma_bwd$degree))
  polm_ma_fwd4cpp = matrix(polm_ma_fwd[,, iseq(2, 1+tmpl$ma_fwd$degree)], 
                           nrow = dim_out, ncol = dim_out*tmpl$ma_fwd$degree)
  
  data_out = matrix(0, 
                    nrow = dim_out, ncol = n_obs)
  
  data_out = get_residuals(data_in = t(data_long), 
                           polm_ar = polm_ar4cpp, polm_ma_bwd = polm_ma_bwd4cpp, polm_ma_fwd = polm_ma_fwd4cpp, 
                           kappa = kappa, k = k)
  
  data_out = data_out[, (tmpl$ar$degree+tmpl$ma_bwd$degree+2):(n_obs-tmpl$ma_fwd$degree)]
  n_valid = dim(data_out)[2]
  
  cov = tcrossprod(data_out)/n_valid
  cov_eigen = eigen(cov)
  
  B = cov_eigen$vectors %*% diag(sqrt(cov_eigen$values)) %*% solve(cov_eigen$vectors)
  
  # Noise Distribution
  distr_par =
    if (tmpl$distr$shock_distr %in% c("gaussian", "laplace")){
      double(0)
    } else if (tmpl$distr$shock_distr %in% c("mixture", "sgt", "skewed_t")) {
      stop("get_init_armamod_whf_random(): This function should not be called for distributions of type *mixture* or *sgt")
    } else {
      stop("get_init_armamod_whf_random(): This should not happen.")
    }
  
  param_noise = c(c(B),
                  c(distr_par))
  
  return(c(param_sys, param_noise))  
}