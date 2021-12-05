#' Function factory for gradient of log-likelihood function
#' 
#' Not used for optimisation but rather obtaining the OPG version of the information matrix.
#' \cr
#' Creates a closure (function plus enclosing environment) which depends only on the deep parameters.
#' In this way, for each template, i.e. integer-valued parameters, and each shock distribution, a log-likelihood function can be created which also contains the data set.
#' It uses similar code as \link{ll_whf_factory} but only exists for the SGT and the skewed t distribution.
#' 
#' @param data_wide Matrix. Data in wide form, i.e. every column contains an observation.
#' @param tmpl List. Containing the mappings from deep to linear parameters, see return value of \link{tmpl_whf_rev}.
#' @param use_cpp Boolean. Default TRUE. If FALSE, then the R version is used instead of the faster Rcpp version.
#'   Used for debugging and checking.
#' @param shock_distr Default \code{sgt}. Other option is \code{skewed_t}
#'
#' @return Closure for particular SVARMA-WHF template \code{tmpl} using the SGT distribution.
#'   Depends only on deep parameters, the dataset is contained in the enclosing environment.
#' 
#' @export
ll_whf_factory_grad_helper = function(data_wide, tmpl, use_cpp = TRUE, shock_distr = "sgt"){
  
  
  # Check inputs and integer-valued parameters ####
  
  # Data
  stopifnot("Input argument *data_wide* must be a matrix." = is.matrix(data_wide),
            "Input argument *data_wide* must contain each observation in a column, 
             and there must be more observations than variables" = dim(data_wide)[2] > dim(data_wide)[1])
  dim_out = dim(data_wide)[1]
  n_obs = dim(data_wide)[2]
  
  # Number of noise parameters in theta
  n_noise_B = dim_out^2
  if (shock_distr == "sgt"){
    n_noise_distr = double(3*dim_out)
  } else if (shock_distr == "skewed_t"){
    n_noise_distr = double(2*dim_out)
  } else {
    stop("ll_whf_factory_grad_helper(): something went wrong with shock_distr argument")
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
    
    # Transform to shocks ####
    B = armamod_whf_obj$B
    eigen_B = eigen(B, only.values = TRUE)$values
    
    data_out = solve(B, data_out)
    
    # Distribution ####
    if (shock_distr == "sgt"){
      distr_params_mat = matrix(theta[(n_par_sys+dim_out^2+1):(n_par_sys+dim_out^2+3*dim_out)],
                                nrow = 3, ncol = dim_out)
      lambda_all = distr_params_mat[1, ]
      p_all = distr_params_mat[2, ]
      q_all = distr_params_mat[3, ]
      
      data_out_tb_cols = purrr::map_dfc(1:dim_out, 
                                        ~ sgt::dsgt(x = data_out[.x, ], 
                                                    mu = 0, sigma = 1,
                                                    mean.cent = TRUE, var.adj = TRUE,
                                                    lambda = lambda_all[.x],
                                                    p = p_all[.x],
                                                    q = q_all[.x],
                                                    log = TRUE) %>% 
                                          dplyr::as_tibble() %>% 
                                          purrr::set_names(nm = paste0("var_", .x)))
    } else if (shock_distr == "skewed_t"){
      distr_params_mat = matrix(theta[(n_par_sys+dim_out^2+1):(n_par_sys+dim_out^2+2*dim_out)],
                                nrow = 2, ncol = dim_out)
      lambda_all = distr_params_mat[1, ]
      q_all = distr_params_mat[2, ]
      
      data_out_tb_cols = purrr::map_dfc(1:dim_out, 
                                        ~ sgt::dsgt(x = data_out[.x, ], 
                                                    mu = 0, sigma = 1,
                                                    mean.cent = TRUE, var.adj = TRUE,
                                                    lambda = lambda_all[.x],
                                                    p = 2,
                                                    q = q_all[.x],
                                                    log = TRUE) %>% 
                                          dplyr::as_tibble() %>% 
                                          purrr::set_names(nm = paste0("var_", .x)))
    } else {
      stop("ll_whf_factory_grad_helper(): specify shock_distr argument correctly")
    }
    
    
    data_out_tb_cols = data_out_tb_cols %>% 
      dplyr::mutate(sum_vars = rowSums(dplyr::across(tidyselect::starts_with("var_")))) %>% 
      tibble::add_column(noise_part = - sum(log(abs(eigen_B)))) %>% 
      dplyr::mutate(result = .data$sum_vars + .data$noise_part)

    return(data_out_tb_cols %>% dplyr::pull(result))
  }
  
  return(ll_whf)
}
