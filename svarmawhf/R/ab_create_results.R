#' Create results
#' 
#' Called by main script \code{scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R}.
#' Calls \link{ll_whf_factory}, \link{get_ic}, \link{replace_noise}.
#' \cr
#' This function is called in a parallelized fashion by the main script \code{scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R} (which in turn is called by the SLURM batch job, see file \code{scriptsP_whf_revision_sgt_bq/b_slurm}.
#' \cr
#' For given initial values and a template, it creates all (optimisation) output which is needed for further analysis, see the return value below.
#'
#' @param theta_init vector of doubles. Initial deep parameters. Obtained from \link{get_init_armamod_whf_random}.
#' @param tmpl SVARMA-WHF template. Obtained from \link{tmpl_whf_rev}.
#' @param params List. Contains optimisation parameters, some of which are taken from the SLURM script which calls this function. 
#'   The SLURM slots are
#'   \itemize{
#'     \item \code{USE_PARALLEL}: Actually set in the main script. Set to FALSE because we use array/jobs \url{https://scicomp.aalto.fi/triton/tut/array/} such that each job only uses one core. This is the most efficient way for embarrassingly parallel problems.
#'     \item \code{N_CORES}: Only relevant if slot \code{USE_PARALLEL} is set to TRUE. In that case, it is retrieved from the SLURM system environment.
#'     \item \code{N_MODS_PER_CORE}: Important parameter. Specifies how many models are estimated by each array-job. 
#'       Given as first parameter \code{$N_MODS_PER_CORE} in SLURM script \code{scriptsP_whf_revision_sgt_bq/b_slurm}.
#'     \item \code{IX_ARRAY_JOB}: Index of array-job. Number of array-jobs is determined from number of rows of dataframe containing all integer-values parameters
#'       Given as second parameter \code{$SLURM_ARRAY_TASK_ID} in SLURM script \code{scriptsP_whf_revision_sgt_bq/b_slurm}.
#'     \item \code{SLURM_JOB_ID}: Internal ID.
#'       Given as third parameter \code{$SLURM_JOB_ID} in SLURM script \code{scriptsP_whf_revision_sgt_bq/b_slurm}.
#'     \item \code{MANUALLY_ASSIGNED_ID}: Manually assigned job ID (same for each array job). Used to keep track of different runs.
#'       Given as fourth parameter \code{$MANUALLY_ASSIGNED_ID} in SLURM script \code{scriptsP_whf_revision_sgt_bq/b_slurm}.
#'   }
#'   The \emph{filename} slots are
#'   \itemize{
#'     \item \code{FILE_NAME_INPUT}: Path to rds file containing the data set. 
#'       Working directory is the location of the main script \code{scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R} called by \code{scriptsP_whf_revision_sgt_bq/b_slurm}.
#'       It is set to \code{"../local_data/g_gmr_bq/data_xts.rds"} in \code{scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R}.
#'     \item \code{PATH_RESULTS_HELPER}: Helper for creating the filenames of the outputs (rds files for each array job). 
#'       Filenames will eventually involve \code{MANUALLY_ASSIGNED_ID} and \code{IX_ARRAY_JOB}.
#'       Set to \code{"../local_data/p_whf/ukko_bq_"} in \code{scriptsP_whf_revision_sgt_bq/aa_Rscript4slurm_bq.R}.
#'   }
#'   The slots for maximal AR and MA orders (for each combination of (p,q) all partial indices are obtained) are:
#'   \itemize{
#'     \item \code{AR_ORDER_MAX}
#'     \item \code{MA_ORDER_MAX}
#'   }
#'   Next, there are some slots regarding optimisation. 
#'   It is advantageous to restart optimisation rather than increase the number of maximal iterations in \link[stats]{optim}.
#'   Also, using different methods (Nelder-Mead vs BFGS) has advantages over just using one method.
#'   For each restart, we perform first a BFGS optimisation, followed by Nelder-Mead.
#'   Only if the value of the likelihood function improves the corresponding optimizing argument is used in the subsequent optimisation.
#'   In order to investigate the convergence behaviour we perform these steps for different shock densities (Gaussian, Laplace, SGT).
#'   \itemize{
#'     \item \code{IT_OPTIM_GAUSS}: Number of restarts for Gaussian density. Set to 3.
#'     \item \code{USE_BFGS_GAUSS}: boolean; whether BFGS method should be used in \link[stats]{optim}. Set to \code{TRUE}.
#'     \item \code{USE_NM_GAUSS}: boolean; whether BFGS method should be used in \link[stats]{optim}. Set to \code{TRUE}.
#'     \item \code{MAXIT_BFGS_GAUSS}: integer. Number of maximal iterations in BFGS method. Set to 80. Default is 100 in \link[stats]{optim}.
#'     \item \code{MAXIT_NM_GAUSS}: integer. Number of maximal iterations in Nelder-Mead method. Set to 1000. Default is 500 in \link[stats]{optim}.
#'     \item \code{IT_OPTIM_LAPLACE}: Same as above but for Laplace density. Set to 3.
#'     \item \code{USE_BFGS_LAPLACE}: Same as above but for Laplace density. Set to \code{TRUE}.
#'     \item \code{USE_NM_LAPLACE}: Same as above but for Laplace density. Set to \code{TRUE}.
#'     \item \code{MAXIT_BFGS_LAPLACE}: Same as above but for Laplace density. Set to 100. 
#'     \item \code{MAXIT_NM_LAPLACE}: Same as above but for Laplace density. Set to 2000.
#'     \item \code{IT_OPTIM_SGT}: Same as above but for SGT densities. Set to 4.
#'     \item \code{USE_BFGS_SGT}: Same as above but for SGT densities. Set to \code{TRUE}.
#'     \item \code{USE_NM_SGT}: Same as above but for SGT densities. Set to \code{TRUE}.
#'     \item \code{MAXIT_BFGS_SGT}: Same as above but for SGT densities. Set to 100.
#'     \item \code{MAXIT_NM_SGT}: Same as above but for SGT densities. Set to 3000.
#'  }
#' @param DATASET matrix. data where each column corresponds to a variable and each row to an observation.
#'
#' @return List with slots \code{params_deep_final}, containing the optimizing deep parameter for the SGT density, \code{value_final}, value of the optimisation function, \code{results_list}, list containing detailed output.
#' The slot \code{results_list} has, in turn, the following slots:
#' \itemize{
#' \item \code{scriptparams}: equal to input argument \code{params}.
#' \item \code{input_integerparams}: (p,q, kappa, k)
#' \item \code{initial}: List with slots \code{theta} (for deep parameters) and \code{value_laplace} (value of laplace density for the initial deep parameters, which is tracked throughout all optimisations including the ones for the other densities).
#' \item \code{gaussian}: List which contains all results for the Gaussian density optimisation.
#'   Number of (unnamed) slots is equal to \code{params$IT_OPTIM_GAUSS}. 
#'   Each of these slots contains a list with slots \code{BFGS} and \code{NM} which are in turn lists that store the output of \link[stats]{optim}. 
#'   The slots for each optimisation method contain the following named slots:
#'   \itemize{
#'     \item \code{convergence}: Integer describing the convergence behaviour of the optimisation in \link[stats]{optim}. 
#'       0 means convergence, 1 means number of maximal iterations reached before convergence (this happens more often for Nelder-Mead than for BFGS), 
#'       2 is a code defined by myself and indicates whether \link[stats]{optim} threw an error, 10 stands for a degenerate NM simplex.
#'     \item \code{duration}: saves the runtime of the optimisation
#'     \item \code{msg}: Error message or success message
#'     \item \code{theta}: deep parameters
#'     \item \code{value}: value of likelihood function of the current kind
#'     \item \code{value_laplace}: value of Laplace likelihood which is tracked throughout all optimisations
#'   }
#' \item \code{ic} list with slots \code{theta} (after rotating the noise parameters using \link[steadyICA]{steadyICA}), and \code{value_laplace}. 
#' \item \code{laplace} and \code{sgt}: Same as slot \code{gaussian} but for Laplace density and SGT densities
#' }
#' 
#' @export
create_results_list = function(theta_init, tmpl, params, DATASET){
 
  # Create some functions (safely) ####
  optim_s = purrr::safely(stats::optim)
  
  ll_fun_gaussian = ll_whf_factory(data_wide = t(DATASET), 
                                   tmpl = tmpl, shock_distr = "gaussian", 
                                   use_cpp = TRUE)
  ll_fun_laplace = ll_whf_factory(data_wide = t(DATASET), 
                                  tmpl = tmpl, shock_distr = "laplace", 
                                  use_cpp = TRUE)
  
  results = list()
  results$script_params = params
  results$input_integerparams = list(p = tmpl$input_orders$ARorder,
                                     q = tmpl$input_orders$MAorder,
                                     kappa = tmpl$input_orders$kappa,
                                     k = tmpl$input_orders$k)
  print(results$input_integerparams)
  
  # Initial values (only storing, obtained elsewhere) ####
  results$initial$theta = theta_init
  results$initial$value_laplace = ll_fun_laplace(results$initial$theta)
  
  val0 = results$initial$value_laplace
  th0 = results$initial$theta
  
  # Gaussian ####
  results$gaussian = vector("list", params$IT_OPTIM_GAUSS)
  for(ix_run in 1:params$IT_OPTIM_GAUSS){
    cat(paste0("Gaussian optimization, iteration number: ", ix_run, "\n"))
    
    if (params$USE_BFGS_GAUSS){
      start_time = Sys.time()
      
      optim_out = optim_s(th0, ll_fun_gaussian, 
                          method = "BFGS", control = list(maxit = params$MAXIT_BFGS_GAUSS))
      
      end_time = Sys.time()
      
      if(is.null(optim_out$result)){
        results$gaussian[[ix_run]]$BFGS$convergence = 2
        results$gaussian[[ix_run]]$BFGS$duration = end_time-start_time
        results$gaussian[[ix_run]]$BFGS$msg = optim_out$error
        results$gaussian[[ix_run]]$BFGS$theta = rep(NA_real_, length(th0))
        results$gaussian[[ix_run]]$BFGS$value = Inf
        results$gaussian[[ix_run]]$BFGS$value_laplace = Inf
      } else {
        results$gaussian[[ix_run]]$BFGS$convergence = optim_out$result$convergence
        results$gaussian[[ix_run]]$BFGS$duration = end_time-start_time
        results$gaussian[[ix_run]]$BFGS$msg = optim_out$result$message
        results$gaussian[[ix_run]]$BFGS$theta = optim_out$result$par
        results$gaussian[[ix_run]]$BFGS$value = optim_out$result$value
        results$gaussian[[ix_run]]$BFGS$value_laplace = ll_fun_laplace(optim_out$result$par)
        
        if(optim_out$result$value < val0){
          th0 <- optim_out$result$par
          val0 <- optim_out$result$value
        }
      }
    }
    
    if (params$USE_NM_GAUSS){
      start_time = Sys.time()
      
      optim_out = optim_s(th0, ll_fun_gaussian, 
                          method = "Nelder-Mead", control = list(maxit = params$MAXIT_NM_GAUSS))
      
      end_time = Sys.time()
      
      if(is.null(optim_out$result)){
        results$gaussian[[ix_run]]$NM$convergence = 2
        results$gaussian[[ix_run]]$NM$duration = end_time-start_time
        results$gaussian[[ix_run]]$NM$msg = optim_out$error
        results$gaussian[[ix_run]]$NM$theta = rep(NA_real_, length(th0))
        results$gaussian[[ix_run]]$NM$value = Inf
        results$gaussian[[ix_run]]$NM$value_laplace = Inf
      } else {
        results$gaussian[[ix_run]]$NM$convergence = optim_out$result$convergence
        results$gaussian[[ix_run]]$NM$duration = end_time-start_time
        results$gaussian[[ix_run]]$NM$msg = optim_out$result$message
        results$gaussian[[ix_run]]$NM$theta = optim_out$result$par
        results$gaussian[[ix_run]]$NM$value = optim_out$result$value
        results$gaussian[[ix_run]]$NM$value_laplace = ll_fun_laplace(optim_out$result$par)
        
        if(optim_out$result$value < val0){
          th0 <- optim_out$result$par
          val0 <- optim_out$result$value
        }
      }
    }
    
  }
  
  rm(ll_fun_gaussian)
  
  # IC ####
  results$ic$theta = replace_noise(th0, 
                                   get_ic(get_residuals_once(params_deep = th0, 
                                                             tmpl = tmpl, 
                                                             DATASET)))
  results$ic$value_laplace = ll_fun_laplace(results$ic$theta)
  
  if(is.numeric(results$ic$theta)){
    th0 = results$ic$theta
    val0 = results$ic$value_laplace
  }
  
  # Laplace ####
  results$laplace = vector("list", params$IT_OPTIM_LAPLACE)
  for(ix_run in 1:params$IT_OPTIM_LAPLACE){
    cat(paste0("Laplace optimization, iteration number: ", ix_run, "\n"))
    
    if (params$USE_BFGS_LAPLACE){
      start_time = Sys.time()
      
      optim_out = optim_s(th0, ll_fun_laplace, 
                          method = "BFGS", control = list(maxit = params$MAXIT_BFGS_LAPLACE))
      
      end_time = Sys.time()
      
      if(is.null(optim_out$result)){
        results$laplace[[ix_run]]$BFGS$convergence = 2
        results$laplace[[ix_run]]$BFGS$duration = end_time-start_time
        results$laplace[[ix_run]]$BFGS$msg = optim_out$error
        results$laplace[[ix_run]]$BFGS$theta = rep(NA_real_, length(th0))
        results$laplace[[ix_run]]$BFGS$value = Inf
        results$laplace[[ix_run]]$BFGS$value_laplace = Inf
      } else {
        results$laplace[[ix_run]]$BFGS$convergence = optim_out$result$convergence
        results$laplace[[ix_run]]$BFGS$duration = end_time-start_time
        results$laplace[[ix_run]]$BFGS$msg = optim_out$result$message
        results$laplace[[ix_run]]$BFGS$theta = optim_out$result$par
        results$laplace[[ix_run]]$BFGS$value = optim_out$result$value
        results$laplace[[ix_run]]$BFGS$value_laplace = ll_fun_laplace(optim_out$result$par)
        
        if(optim_out$result$value < val0){
          th0 <- optim_out$result$par
          val0 <- optim_out$result$value
        }
      }
    }
    
    if (params$USE_NM_LAPLACE){
      start_time = Sys.time()
      
      optim_out = optim_s(th0, ll_fun_laplace, 
                          method = "Nelder-Mead", control = list(maxit = params$MAXIT_NM_LAPLACE))
      
      end_time = Sys.time()
      
      if(is.null(optim_out$result)){
        results$laplace[[ix_run]]$NM$convergence = 2
        results$laplace[[ix_run]]$NM$duration = end_time-start_time
        results$laplace[[ix_run]]$NM$msg = optim_out$error
        results$laplace[[ix_run]]$NM$theta = rep(NA_real_, length(th0))
        results$laplace[[ix_run]]$NM$value = Inf
        results$laplace[[ix_run]]$NM$value_laplace = Inf
      } else {
        results$laplace[[ix_run]]$NM$convergence = optim_out$result$convergence
        results$laplace[[ix_run]]$NM$duration = end_time-start_time
        results$laplace[[ix_run]]$NM$msg = optim_out$result$message
        results$laplace[[ix_run]]$NM$theta = optim_out$result$par
        results$laplace[[ix_run]]$NM$value = optim_out$result$value
        results$laplace[[ix_run]]$NM$value_laplace = ll_fun_laplace(optim_out$result$par)
        
        if(optim_out$result$value < val0){
          th0 <- optim_out$result$par
          val0 <- optim_out$result$value
        }
      }
    }
    
  }
  
  # Initial SGT parameter estimation ####
  res4sgt = get_residuals_once(th0, tmpl, DATASET)
  
  # safely version of fitdist: If error, return starting values
  
  # apply fct for creating distribution parameters
  
  fitdist_s = purrr::safely(fitdistrplus::fitdist, otherwise = list(estimate = c(0,1.5,20)))
  
  fitdist_obj = purrr::map(1:params$DIM_OUT, ~ suppressWarnings(fitdist_s(data = res4sgt[,.x],
                                                                   distr = "sgt", 
                                                                   start = list(lambda = 0, p = 1.5, q = 20),
                                                                   fix.arg = list(mu = 0, sigma = 1), 
                                                                   discrete = FALSE)$result$estimate))
  
  distr_init_params = matrix(NA, nrow = 3, ncol = params$DIM_OUT)
  for (ix_comp in 1:params$DIM_OUT){
    distr_init_params[,ix_comp] = fitdist_obj[[ix_comp]]
  }
  
  rm(res4sgt)
  n_nondistr_hlp = length(th0)
  n_distr_hlp = length(c(distr_init_params))
  th0 = c(th0, c(distr_init_params))
  
  # SGT log-likelihood ####
  
  # Update template for distribution
  tmpl$distr = list(n_par = params$DIM_OUT * 3,
                    shock_distr = "sgt")
  
  ll_fun_sgt = ll_whf_factory(data_wide = t(DATASET), 
                              tmpl = tmpl, shock_distr = "sgt", 
                              use_cpp = TRUE)
  val0 = ll_fun_sgt(th0)
  
  # SGT ####
  results$sgt = vector("list", params$IT_OPTIM_SGT)
  for(ix_run in 1:params$IT_OPTIM_SGT){
    cat(paste0("SGT optimization, iteration number: ", ix_run, "\n"))
    
    if (params$USE_BFGS_SGT){
      start_time = Sys.time()
      
      optim_out = optim_s(th0, ll_fun_sgt, 
                          method = "L-BFGS-B", 
                          lower = c(rep(-Inf, n_nondistr_hlp), rep(c(-1, 1.5, 2), params$DIM_OUT)),
                          upper = c(rep(Inf, n_nondistr_hlp), rep(c(1, Inf, Inf), params$DIM_OUT)),
                          control = list(maxit = params$MAXIT_BFGS_SGT))
      
      end_time = Sys.time()
      
      if(is.null(optim_out$result)){
        results$sgt[[ix_run]]$BFGS$convergence = 2
        results$sgt[[ix_run]]$BFGS$duration = end_time-start_time
        results$sgt[[ix_run]]$BFGS$msg = optim_out$error
        results$sgt[[ix_run]]$BFGS$theta = rep(NA_real_, length(th0))
        results$sgt[[ix_run]]$BFGS$value = Inf
        results$sgt[[ix_run]]$BFGS$value_laplace = Inf
      } else {
        results$sgt[[ix_run]]$BFGS$convergence = optim_out$result$convergence
        results$sgt[[ix_run]]$BFGS$duration = end_time-start_time
        results$sgt[[ix_run]]$BFGS$msg = optim_out$result$message
        results$sgt[[ix_run]]$BFGS$theta = optim_out$result$par
        results$sgt[[ix_run]]$BFGS$value = optim_out$result$value
        results$sgt[[ix_run]]$BFGS$value_laplace = ll_fun_laplace(optim_out$result$par[1:(length(th0)-tmpl$distr$n_par)])
        
        if(optim_out$result$value < val0){
          th0 <- optim_out$result$par
          val0 <- optim_out$result$value
        }
      }
    }
    
    if (params$USE_NM_SGT){
      start_time = Sys.time()
      
      optim_out = optim_s(th0, ll_fun_sgt, 
                          method = "Nelder-Mead", control = list(maxit = params$MAXIT_NM_SGT))
      
      end_time = Sys.time()
      
      if(is.null(optim_out$result)){
        results$sgt[[ix_run]]$NM$convergence = 2
        results$sgt[[ix_run]]$NM$duration = end_time-start_time
        results$sgt[[ix_run]]$NM$msg = optim_out$error
        results$sgt[[ix_run]]$NM$theta = rep(NA_real_, length(th0))
        results$sgt[[ix_run]]$NM$value = Inf
        results$sgt[[ix_run]]$NM$value_laplace = Inf
      } else {
        results$sgt[[ix_run]]$NM$convergence = optim_out$result$convergence
        results$sgt[[ix_run]]$NM$duration = end_time-start_time
        results$sgt[[ix_run]]$NM$msg = optim_out$result$message
        results$sgt[[ix_run]]$NM$theta = optim_out$result$par
        results$sgt[[ix_run]]$NM$value = optim_out$result$value
        results$sgt[[ix_run]]$NM$value_laplace = ll_fun_laplace(optim_out$result$par[1:(length(th0)-tmpl$distr$n_par)])
        
        if(optim_out$result$value < val0){
          th0 <- optim_out$result$par
          val0 <- optim_out$result$value
        }
      }
    }
    
  }
  
  rm(ll_fun_laplace) # only removed here because it also tracks the likelihood values in SGT optimization
  
  return(list(params_deep_final = th0,
              value_final = val0,
              results_list = results))
}
