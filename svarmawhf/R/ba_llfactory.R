# Used for Creating Initial Values ----

#' Extract deep parameters from SARMA-WHF model given a template
#' 
#' Called by \link{get_init_armamod_whf_random}.
#' \cr
#' Uses the affine mappings in \code{tmpl} to reconstruct the deep parameters from the linear parameters in \code{armamod_whf_obj} with \link[stats]{lsfit}.
#'
#' @param armamod_whf_obj List of 3D arrays \code{polm_ar}, \code{polm_ma_bwd}, \code{polm_ma_fwd}, and the static shock transmission matrix \code{B}, 
#'   as well as info regarding the distribution in slot \code{distr} Obtained from \link{fill_tmpl_whf_rev}.
#' @param tmpl List. Affine mappings etc specifiying a SVARMA-WHF template. See return value of \link{tmpl_whf_rev}.
#' @param tol Tolerance level. "Residuals" when extracting deep parameters must be smaller than this.
#'   Otherwise an error is thrown.
#'   Default equal to \code{sqrt(.Machine$double.eps)}
#'
#' @return Vector of doubles containing deep parameters.
#' @export
extract_deep_rev = function(armamod_whf_obj, tmpl, tol = sqrt(.Machine$double.eps)){
  
  mod_whf_vec = c(c(armamod_whf_obj$polm_ar),
                  c(armamod_whf_obj$polm_ma_bwd),
                  c(armamod_whf_obj$polm_ma_fwd))
  
  h = c(c(tmpl$ar$h),
        c(tmpl$ma_bwd$h),
        c(tmpl$ma_fwd$h))
  
  H = bdiag(tmpl$ar$H,
            tmpl$ma_bwd$H,
            tmpl$ma_fwd$H)
  
  # mod_whf_vec = h + H vec_deep
  out = stats::lsfit(H, mod_whf_vec-h, intercept = FALSE)
  
  stopifnot("extract_deep_rev(): Residuals should be zero when deep parameters are extracted using the ARMA-WHF template." = max(abs(out$residuals)) < tol)
  
  return(unname(out$coef))
}


#' Get Residuals
#' 
#' Calculates the residuals by using \link{toepl_inv} (MA calculations) and \link{toepl_fwd} (AR calculations).
#' The inputs and the way this function is written mimics the Rcpp function which is actually used in the calculations because of speed.
#' In particular, the AR, MA bwd, MA fwd polynomial matrices are not in the usual 3D array format. 
#' Moreover, they differ in whether the zero-lag coefficient matrix is included and in which direction they are ordered (this is due to computational reasons).
#' 
#' @param data_in matrix. Observed data in wide format, i.e. observations in columns and variables in rows.
#'   This is done for consistency with the Rcpp version \link{get_residuals}
#' @param polm_ar matrix. AR polynomial in wide matrix format (not a 3D array). \eqn{( I, a_1, \ldots, a_p )}{(I, a[1], ... , a[p])}
#'     of dimension \eqn{(n \times (p+1) n)}{(n x (p+1)n)}.
#'     It is used in the call to \link{toepl_fwd}.
#' @param polm_ma_bwd matrix. MA backward polynomial in wide matrix format (not a 3D array). 
#'   It includes the zero lag coefficient but is in reverse order: \eqn{p_{q-kappa}, \ldots, p_1, p_0}.
#'   The matrix \eqn{p_0} is not necessarily the identity matrix and this needs to be taken into account before discarding it such that 
#'   \link{toepl_inv} can be called.
#' @param polm_ma_fwd matrix. MA forward polynomial in wide matrix format (not a 3D array). Zero lag coefficient is \strong{not} included.
#'   Note that for usual MA matrix polynomials the order of the coefficients would need to be inverted to call \link{toepl_inv}.
#'   However, because of the forward nature of these calculations, the matrices are in the order 
#'   \eqn{( f_1, \ldots, f_{kappa} )}{(f[1], ... , f[kappa])} and of dimension \eqn{(n \times kappa n)}{(n x kappa n)}.
#' @param kappa,k partial index
#' 
#' @return Residuals in wide format (every row corresponds to a variable, every column to an observation)
#'
#' @export
get_residuals_R = function(data_in, 
                           polm_ar, polm_ma_bwd, polm_ma_fwd, 
                           kappa, k){
  
  # Retrieve integer valued parameters
  dim_out = dim(data_in)[1]
  n_obs = dim(data_in)[2]
  
  ARorder = dim(polm_ar)[2]/dim_out - 1
  MAorder_bwd = dim(polm_ma_bwd)[2]/dim_out - 1
  MAorder_fwd = dim(polm_ma_fwd)[2]/dim_out
  
  # AR calculations w_t = a(z) y_t
  data_ar = toepl_fwd(polm_ar, data_in)
  
  # MA calculations: Solve w_t = p(z) u_t for u_t (and take non-identity p0 into account)
  p0_inv = polm_ma_bwd[, (MAorder_bwd*dim_out + 1):(MAorder_bwd*dim_out + dim_out), drop = FALSE] 
  polm_ma_bwd = p0_inv %*% polm_ma_bwd[, 1:(MAorder_bwd*dim_out), drop = FALSE]
  data_ar = p0_inv %*% data_ar
  
  data_ma_bwd = toepl_inv(polm_ma_bwd, data_ar)
  
  # Shift with k of partial index (note that kappa has no effect on Toeplitz calculations)
  if (k > 0){
    data_ma_bwd[(k+1):dim_out, 2:n_obs] = data_ma_bwd[(k+1):dim_out, 1:(n_obs-1), drop = FALSE]
    data_ma_bwd[,1] = 0
  }
  
  # MA forward calculations: Solve u_t = f(z) v_t for v_t (BUT IN REVERSE TIME ORDER)
  data_ma_fwd = toepl_inv(polm_ma_fwd, data_ma_bwd[, n_obs:1, drop = FALSE])
  
  # Reverse time order back to normal again
  data_ma_fwd = data_ma_fwd[, n_obs:1, drop = FALSE]
  
  return(data_ma_fwd)
}


# Used for Creating Results ----

#' Function factory for negative log-likelihood minimization
#' 
#' Called by \link{create_results_list}.
#' Calls \link{fill_tmpl_whf_rev}, \link{get_residuals}.
#' \cr
#' Creates a closure (function plus enclosing environment) which depends only on the deep parameters.
#' In this way, for each template, i.e. integer-valued parameters, and each shock distribution, a log-likelihood function can be created which also contains the data set.
#' This is convenient for optimisation via \link[stats]{optim}.
#' \cr
#' In the main script \code{scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R}, for each tuple of integer-valued parameters \eqn{(p,q,\kappa,k)} a log-likelihood function for each shock distribution \code{shock_distr} is created.
#'
#' @param data_wide Matrix. Data in wide form, i.e. every column contains an observation.
#' @param tmpl List. Containing the mappings from deep to linear parameters, see return value of \link{tmpl_whf_rev}.
#' @param shock_distr Character. One of \code{laplace}, \code{mixture}, \code{gaussian}, \code{sgt}, \code{skewed_t}, \code{skewed_ged}
#' @param noise_only boolean. Default set to \code{FALSE}. If \code{TRUE}, then in the likelihood function to be created only the noise parameters are optimized (it depends only on the noise parameters in this case).
#' @param use_cpp Boolean. Default TRUE. If FALSE, then the R version is used instead of the faster Rcpp version.
#'   Used for debugging and checking.
#' @param info Boolean. Default set to \code{FALSE}. 
#'   Only used for debugging in which case tibble_info and call_counter must be in the global environment.
#'
#' @return Closure for particular SVARMA-WHF template \code{tmpl} and shock distribution \code{shock_distr} which depends only on deep parameters.
#' 
#' @export
#' 
#' @importFrom LaplacesDemon dnormm
#' @importFrom LaplacesDemon dlaplace
#' @importFrom purrr map_dbl
ll_whf_factory = function(data_wide, tmpl, shock_distr = c("laplace", "mixture", "gaussian", "sgt", "skewed_t", "skewed_ged"), noise_only = FALSE, use_cpp = TRUE, info = FALSE){
  
  # promises
  force(data_wide)
  force(tmpl)
  shock_distr = match.arg(shock_distr)
  force(shock_distr)
  
  # Check inputs and integer-valued parameters ####
  
  # Data
  stopifnot("Input argument *data_wide* must be a matrix." = is.matrix(data_wide),
            "Input argument *data_wide* must contain each observation in a column, 
             and there must be more observations than variables" = dim(data_wide)[2] > dim(data_wide)[1])
  dim_out = dim(data_wide)[1]
  n_obs = dim(data_wide)[2]
  
  # Number of noise parameters in theta
  n_noise_B = dim_out^2
  n_noise_distr = if (shock_distr %in% c("laplace", "gaussian")){
    double(0L)
  } else if (shock_distr %in% c("mixture", "sgt")){
    double(3*dim_out)
  } else if (shock_distr == "skewed_t"){
    double(2*dim_out)
  } else {
    stop("ll_whf_factory(): Choose an available *shock_distribution*")
  }
  n_noise = n_noise_B + n_noise_distr
  n_par_sys = tmpl$ar$n_par + tmpl$ma_bwd$n_par + tmpl$ma_fwd$n_par
  

  ll_whf = function(theta){
    
    # Construct system parameters from template ####
    armamod_whf_obj = fill_tmpl_whf_rev(theta, tmpl)
    
    polm_ar = armamod_whf_obj$polm_ar
    polm_ma_bwd = armamod_whf_obj$polm_ma_bwd
    polm_ma_fwd = armamod_whf_obj$polm_ma_fwd
    
    # Check eigenvalues ####
    if (tmpl$ar$degree > 0){
      eigen_ar = eigen(companion_matrix(polm_ar), only.values = TRUE)$values
    } else {
      eigen_ar = 0
    }
    
    if (tmpl$ma_bwd$degree > 0){
      eigen_ma_bwd = eigen(companion_matrix(polm_ma_bwd), only.values = TRUE)$values
    } else {
      eigen_ma_bwd = 0
    }
    
    if (tmpl$ma_fwd$degree > 0){
      eigen_ma_fwd = eigen(companion_matrix(polm_ma_fwd), only.values = TRUE)$values
    } else {
      eigen_ma_fwd = 0
    }
    
    if (any(c(abs(eigen_ar), abs(eigen_ma_bwd), abs(eigen_ma_fwd)) >= 1)) {
      return(1e25)
    }
    
    # Input preparation for RcppArmadillo function get_residuals() ####
    polm_ar4cpp = matrix(polm_ar, 
                         nrow = dim_out, ncol = dim_out*(1 + tmpl$ar$degree))
    polm_ma_bwd4cpp = matrix(polm_ma_bwd[,, rev(iseq(1, 1+tmpl$ma_bwd$degree))], 
                             nrow = dim_out, ncol = dim_out*(1+tmpl$ma_bwd$degree))
    polm_ma_fwd4cpp = matrix(polm_ma_fwd[,, iseq(2, 1+tmpl$ma_fwd$degree)], 
                             nrow = dim_out, ncol = dim_out*tmpl$ma_fwd$degree)
    
    if (use_cpp){
      # data_out = matrix(0, 
      #                   nrow = dim_out, ncol = n_obs)
      
      data_out = get_residuals(data_in = data_wide, 
                               polm_ar = polm_ar4cpp, polm_ma_bwd = polm_ma_bwd4cpp, polm_ma_fwd = polm_ma_fwd4cpp, 
                               kappa = tmpl$input_orders$kappa, k = tmpl$input_orders$k) #, data_out = data_out)
    } else {
      data_out = get_residuals_R(data_in = data_wide, 
                                 polm_ar = polm_ar4cpp, polm_ma_bwd = polm_ma_bwd4cpp, polm_ma_fwd = polm_ma_fwd4cpp, 
                                 kappa = tmpl$input_orders$kappa, k = tmpl$input_orders$k)
    }
    
    data_out = data_out[, (tmpl$ar$degree+tmpl$ma_bwd$degree+2):(n_obs-tmpl$ma_fwd$degree)]
    n_valid = dim(data_out)[2]
    
    B = armamod_whf_obj$B
    eigen_B = eigen(B, only.values = TRUE)$values
    
    data_out = solve(B, data_out)

    # Density specification ####
    if (shock_distr == "mixture"){
      distr_params_mat = matrix(theta[(n_par_sys+dim_out^2+1):(n_par_sys+dim_out^2+3*dim_out)],
                                nrow = 3, ncol = dim_out)
      p_all = distr_params_mat[1, ]
      r_all = distr_params_mat[2, ]
      phi_all = distr_params_mat[3, ]
      
      mu1_all = sqrt((1-p_all)/p_all) * r_all * sin(phi_all)
      mu2_all = -sqrt((1-p_all)/p_all) * r_all * sin(phi_all)
      
      sd1_all = sqrt(1/p_all) * r_all * cos(phi_all)
      sd2_all = sqrt(1/(1-p_all) * (1-r_all))
      
      # parametrise in terms of (p, r, angle)
      ll_part_density = 
        map_dbl(1:dim_out, 
                ~mean(LaplacesDemon::dnormm(x = data_out[.x, ], 
                                            p = c(p_all[.x], 1-p_all[.x]), 
                                            mu = c(mu1_all[.x], mu2_all[.x]),
                                            sigma = c(sd1_all[.x], sd2_all[.x]),
                                            log = TRUE))) %>% 
        sum()  
    } else if (shock_distr == "laplace"){
      ll_part_density = 
        dim_out * mean(LaplacesDemon::dlaplace(data_out, 
                                               location = 0, 
                                               scale = 1/sqrt(2),
                                               log = TRUE))
    } else if (shock_distr == "skewed_t"){
      distr_params_mat = matrix(theta[(n_par_sys+dim_out^2+1):(n_par_sys+dim_out^2+2*dim_out)],
                                nrow = 2, ncol = dim_out)
      lambda_all = distr_params_mat[1, ]
      q_all = distr_params_mat[2, ]
      
      ll_part_density = 
        map_dbl(1:dim_out, 
                ~mean(sgt::dsgt(x = data_out[.x, ], 
                                mu = 0, sigma = 1,
                                mean.cent = TRUE, var.adj = TRUE,
                                lambda = lambda_all[.x],
                                p = 2,
                                q = q_all[.x],
                                log = TRUE))) %>% 
        sum()
    } else if (shock_distr == "skewed_ged"){
      distr_params_mat = matrix(theta[(n_par_sys+dim_out^2+1):(n_par_sys+dim_out^2+2*dim_out)],
                                nrow = 2, ncol = dim_out)
      lambda_all = distr_params_mat[1, ]
      p_all = distr_params_mat[2, ]
      
      ll_part_density = 
        map_dbl(1:dim_out, 
                ~mean(sgt::dsgt(x = data_out[.x, ], 
                                mu = 0, sigma = 1,
                                mean.cent = TRUE, var.adj = TRUE,
                                lambda = lambda_all[.x],
                                p = q_all[.x],
                                q = Inf,
                                log = TRUE))) %>% 
        sum()
    } else if (shock_distr == "sgt"){
      distr_params_mat = matrix(theta[(n_par_sys+dim_out^2+1):(n_par_sys+dim_out^2+3*dim_out)],
                                nrow = 3, ncol = dim_out)
      lambda_all = distr_params_mat[1, ]
      p_all = distr_params_mat[2, ]
      q_all = distr_params_mat[3, ]
      
      ll_part_density = 
        map_dbl(1:dim_out, 
                ~mean(sgt::dsgt(x = data_out[.x, ], 
                                mu = 0, sigma = 1,
                                mean.cent = TRUE, var.adj = TRUE,
                                lambda = lambda_all[.x],
                                p = p_all[.x],
                                q = q_all[.x],
                                log = TRUE))) %>% 
        sum()
     } else {
      ll_part_density = 
        -1/2 * 1/n_valid * crossprod(c(data_out))
    }
    
    if (info){
      # More output for tracing optimization procedure (ugly, but asfaik the only way to do this in R) ####
      tibble_info <<- tibble_info %>% 
        tibble::add_row(call_counter = call_counter,
                        eigen_ar_info = list(eigen_ar),
                        eigen_ma_bwd_info = list(eigen_ma_bwd),
                        eigen_ma_fwd_info = list(eigen_ma_bwd),
                        eigen_B_info = list(eigen_B),
                        ll_info = list(c(-ll_part_density, sum(log(abs(eigen_B))), -ll_part_density + sum(log(abs(eigen_B))))))
      
      call_counter <<- call_counter+1  
    }
    
    # Minus the log-likelihood because we will minimize -1/n_valid loglik(theta)
    return(-ll_part_density + sum(log(abs(eigen_B))))
  }
  
  return(ll_whf)
}

# _Helper for getting Residuals in create_results_list() ----

#' Residuals from deep parameters
#' 
#' Called in \link{create_results_list}.
#' 
#' Different input arguments than \link{get_residuals} which it calls after input transformation.
#' Used when calculating the residuals for the Laplace (after call to \link{get_ic}) and SGT version for the first time. 
#'
#' @param params_deep vector of doubles. System and noise parameters which fill the SVARMA-WHF template.
#' @param tmpl List. Affine mappings etc specifiying a SVARMA-WHF template. See return value of \link{tmpl_whf_rev}.
#' @param data_long matrix. data in long format, i.e. every variable corresponds to a column
#'
#' @return Residuals (long format)
#' @export
get_residuals_once = function(params_deep, tmpl, data_long){
  
  dim_out = dim(data_long)[2]
  
  armamodwhf = fill_tmpl_whf_rev(params_deep, tmpl)
  
  polm_ar4cpp     = matrix(armamodwhf$polm_ar, 
                           nrow = dim_out, ncol = dim_out*(1 + tmpl$ar$degree))
  polm_ma_bwd4cpp = matrix(armamodwhf$polm_ma_bwd[,, rev(iseq(1, 1+tmpl$ma_bwd$degree))], 
                           nrow = dim_out, ncol = dim_out*(1+tmpl$ma_bwd$degree))
  polm_ma_fwd4cpp = matrix(armamodwhf$polm_ma_fwd[,, iseq(2, 1+tmpl$ma_fwd$degree)], 
                           nrow = dim_out, ncol = dim_out*tmpl$ma_fwd$degree)
  
  return(get_residuals(data_in = t(data_long), 
                       polm_ar = polm_ar4cpp, polm_ma_bwd = polm_ma_bwd4cpp, polm_ma_fwd = polm_ma_fwd4cpp, 
                       kappa = tmpl$input_orders$kappa, k = tmpl$input_orders$k) %>% t())
}

# _Initial Rotation of B ----

#' Rotate inputs to Independence
#' 
#' Called in \link{create_results_list}.
#' \cr
#' Uses \link[steadyICA]{steadyICA} in order to obtain independent components.
#' This function is called after the initial Gaussian estimation and should serve as initial estimate for the noise parameters.
#' Afterwards, full ML with Laplace density and SGT densities is performed.
#'
#' @param res Residuals in long form
#'
#' @return Static shock transmission matrix \eqn{B = C M} where \eqn{C} is obtained from 
#'   Cholesky decomposition of \code{1/n_obs t(res) res = C t(C)} and 
#'   \code{M} is the mixing matrix obtained from the output of \link{steadyICA} applied to \eqn{res solve(t(C))}.
#' @export
get_ic = function(res){
  
  res = scale(res, center = TRUE, scale = FALSE)
  
  n_obs = dim(res)[1]
  cov = crossprod(res)/(n_obs-1)
  chol_factor = t(chol(cov))
  
  res_white = t(solve(chol_factor, t(res)))
  
  out_dc = suppressMessages(steadyICA::steadyICA(res_white))
  
  mixing_mat = out_dc$M
  
  return(chol_factor %*% mixing_mat)
}

#' Replace noise parameters
#'
#' Called in \link{create_results_list}.
#' \cr
#' Replaces unidentified estimate obtained from Gaussian estimation with independent components obtained from call to \link{get_ic}.
#' \cr
#' Note that this function makes only sense when shock distribution is Gaussian or Laplace
#' 
#' @param params vector of doubles. Deep system and noise parameters. Its noise parameters will be replaced by \code{params_new_noise}.
#' @param params_new_noise vector of doubles. Noise parameters as obtained from \link{get_ic}.
#'
#' @return Vector of deep parameters \code{params} whose noise-part is updated with rotation such that the shock components are independent.
#' 
#' @export
replace_noise = function(params, params_new_noise){
  n = length(params)
  m = length(c(params_new_noise))
  
  # Only makes sense like this for Gaussian and Laplace density!
  params[(n-m+1):n] = c(params_new_noise)
  return(params)
}



