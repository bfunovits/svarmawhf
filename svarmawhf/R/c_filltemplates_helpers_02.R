# Filling Templates ----

#' Fill WHF template with deep parameters
#' 
#' Called in \link{ll_whf_factory}, \link{armamod_whf}.
#' 
#' Use affine mappings in SVARMA-WHF template to fill free parameters in AR, MA bwd, etc objects (3D arrays) with deep parameters.
#'
#' @param theta Vector of doubles. Deep parameters to be filled in the SVARMA-WHF template \code{tmpl}.
#' @param tmpl SVARMA-WHF template, i.e. list with slots containing affine mappings between deep and linear parameters, see return value of \link{tmpl_whf_rev}.
#'
#' @return List of 3D arrays \code{polm_ar}, \code{polm_ma_bwd}, \code{polm_ma_fwd}, 
#'   and the static shock transmission matrix \code{B}, 
#'   as well as info regarding the distribution in slot \code{distr}
#'   
#' @export
fill_tmpl_whf_rev = function(theta, tmpl){
  
  # Integer-valued parameters
  dim_out = tmpl$input_orders$dim_out
  
  # AR params: Apply affine mapping, then create array by adding a dimension attribute
  polm_ar = tmpl$ar$h + tmpl$ar$H %*% theta[iseq(1, tmpl$ar$n_par)]
  dim(polm_ar) = c(dim_out, dim_out, 1 + tmpl$ar$degree)
  
  # MA bwd params
  polm_ma_bwd = tmpl$ma_bwd$h + tmpl$ma_bwd$H %*% theta[iseq(1+tmpl$ar$n_par, 
                                                             tmpl$ma_bwd$n_par+tmpl$ar$n_par)]
  dim(polm_ma_bwd) = c(dim_out, dim_out, 1 + tmpl$ma_bwd$degree)
  
  # MA fwd params 
  polm_ma_fwd = tmpl$ma_fwd$h + tmpl$ma_fwd$H %*% theta[iseq(1+tmpl$ma_bwd$n_par+tmpl$ar$n_par, 
                                                             tmpl$ma_fwd$n_par+tmpl$ma_bwd$n_par+tmpl$ar$n_par)]
  dim(polm_ma_fwd) = c(dim_out, dim_out, 1 + tmpl$ma_fwd$degree)
  
  # B params
  B = theta[iseq(1 + tmpl$ma_fwd$n_par+tmpl$ma_bwd$n_par+tmpl$ar$n_par, 
                 tmpl$B$n_par + tmpl$ma_fwd$n_par+tmpl$ma_bwd$n_par+tmpl$ar$n_par)]
  dim(B) = c(dim_out, dim_out)
  
  # Distribution
  distr = theta[iseq(tmpl$ma_fwd$n_par+tmpl$ma_bwd$n_par+tmpl$ar$n_par + tmpl$B$n_par + 1,
                     tmpl$ma_fwd$n_par+tmpl$ma_bwd$n_par+tmpl$ar$n_par + tmpl$B$n_par + tmpl$distr$n_par)]
  
  return(list(polm_ar = polm_ar,
              polm_ma_bwd = polm_ma_bwd,
              polm_ma_fwd = polm_ma_fwd,
              B = B,
              distr = distr))
}


#' Generate ARMA Model as (a(z), b(z))
#' 
#' Called when generating IRFs.
#' Calls \link{fill_tmpl_whf_rev}.
#' \cr
#' Takes the deep parameters \code{theta} as input and returns an ARMA model (not in normalised canonical WHF form) 
#' by first filling the ARMA-WHF template \code{tmpl} and subsequent transformation.
#' 
#' @param theta vector of doubles; deep parameters
#' @param tmpl template as obtained by \link{tmpl_whf_rev}, i.e. a list containing slots \code{polm_ar}, \code{polm_ma_bwd}, \code{polm_ma_fwd}, \code{B}.
#' 
#' @return List with slots \code{polm_ar}, \code{polm_ma}, (both \link{polm} objects), and \code{B} (a matrix).
#' @export
armamod_whf = function(theta, tmpl){
  
  # Fill SVARMA-WHF template
  armamodwhf = fill_tmpl_whf_rev(theta, tmpl)
  
  # Transform to polm objects and "usual" MA parameters
  a = armamodwhf$polm_ar %>% polm()
  b1 = armamodwhf$polm_ma_bwd %>% polm()
  b2 = armamodwhf$polm_ma_fwd %>% f2g(c(tmpl$input_orders$kappa, tmpl$input_orders$k))
  B = armamodwhf$B
  
  return(list(polm_ar = a,
              polm_ma = b1 %r% b2,
              B = B))
}

# _Helpers for Filling Templates ----

#' Transform forward polynomial of WHF to polynomial in z
#'
#' Called in \link{armamod_whf}.
#' \cr
#' Given the partial indices, parametrised by \eqn{\kappa} and \code{k},
#' we transform \eqn{f(z)} (a polynomial in \eqn{\frac{1}{z}}) to \eqn{g(z) = s(z)f(z)} and vice versa,
#' see the paper for notation.
#'
#' @param f Array object which contains f(z).
#'   There is no class for polynomials in 1/z so far.
#' @param g \link{polm} object which contains \eqn{g(z) = s(z)f(z)}.
#' @param kappa_k Two dimensional vector of integers parametrising the partial indices of the WHF
#'
#' @return \link{polm} object \eqn{g(z)} or array object for \eqn{f(z)}.
#' @export
f2g = function(f, kappa_k){
  # Checking and integer valued parameters
  f = unclass(f)
  g = f
  
  stopifnot(is.array(f))
  n = dim(f)[1]
  
  kappa = kappa_k[1]
  k = kappa_k[2]
  
  # Shift and reverse coefficients
  if(k==0){
    g = f[,,(kappa+1):1 , drop = FALSE]
  } else {
    g[1:k, ,] = f[1:k, , (kappa+2):1]
    g[(k+1):n, , 1:(kappa+1)] = f[(k+1):n, , (kappa+1):1]
  }
  
  return(polm(g))
}

#' @rdname f2g
#' @export
g2f = function(g, kappa_k){
  f = f2g(g, kappa_k)
  return(unclass(f))
}



# IRF Helpers ----

#' Generate IRF
#' 
#' Called in analysis when comparing IRFs.
#' Calls \link{armamod_whf}.
#' 
#' @param theta,tmpl Used in call to \link{fill_tmpl_whf_rev}.
#' @param n_lags Number of lags for IRF
#'
#' @return \code{pseries} object
#' @export
irf_whf = function(theta, tmpl, n_lags){
  
  armamod_whf = armamod_whf(theta, tmpl)
  
  return(pseries(lmfd(a = armamod_whf$polm_ar, 
                      b = armamod_whf$polm_ma), lag.max = n_lags) %r% armamod_whf$B)
}

#' Global identifiability by choosing a canonical representative
#'
#' Obtaining the canonical representative as described in the identification scheme in Lanne, Meitz, Saikkonen (2017).
#'
#' @param mat A matrix with ones on the diagonal (as obtained from the maximum likelihood procedure).
#' @param sig Vector containing the standard deviations of the i.i.d. shocks  (as obtained from the ML procedure).
#'
#' @return List with items \code{mat} and \code{sig}, containing the permuted and scaled versions of the inputs such that the identification scheme is satisfied.
#'   The diagonal elements of \code{mat} are again equal to one.
#' 
#' @export
#'
#' @examples
#' mm <- matrix(c(1,3,2,1), 2)
#' ss <- c(1,1)
#'
#' ll <- unique_p_s(mm, ss)
#' (m = ll[["mat"]])
#' (s = ll[["sig"]])
#' m %*% diag(s) %*% diag(s) %*% t(m)
#' mm %*% diag(ss) %*% diag(ss) %*% t(mm)
#' 
#' set.seed(1)
#' n = 3
#' mat = matrix(rnorm(n^2), n ,n)
#' unique_p_s(mat)
unique_p_s <- function(mat, sig = NULL){
  
  # Integer-valued params
  dim_out <- ncol(mat)
  
  # Scale columns by their length
  len_col <- apply(mat, MARGIN = 2, function(x) norm(x, "2"))
  scale_1 <- solve(diag(len_col))
  mat <- mat %*% scale_1
  
  # Permute columns such that the diag-element is (weakly) larger than all elements to its right
  perm_mats <- vector("list", dim_out-1)
  for (ix_col in 1:(dim_out-1)){
    tmp_order <- mat[ix_col, ix_col:dim_out, drop = FALSE] %>% abs() %>% order(decreasing = TRUE)
    if (ix_col == 1){
      perm <- tmp_order
    } else{
      perm <- c(1:(ix_col-1), (ix_col-1) + tmp_order )
    }
    tmp_perm <- diag(dim_out)[,perm]
    perm_mats[[ix_col]] <- tmp_perm
    mat <- mat %*% tmp_perm
  }
  
  perm <- Reduce(function(x,y) x %*% y, perm_mats)
  
  # Scale matrix s.t. diagonal elements are equal to one
  scale_2 <- solve(diag(diag(mat)))
  mat <- mat %*% scale_2
  
  # Applying the transformation inv(scale_2) * inv(perm) * inv(scale_1) from left to the diagonal matrix *sig*,
  # we would result in a non-diagonal matrix
  # Note that we can always solve P D = \tilde{D} P such that \tilde{D} is a diagonal matrix
  trafo <- scale_1 %*% perm %*% scale_2
  if (is.null(sig)){
    sig <- trafo[trafo != 0]^(-1)
  } else {
    sig <- sig*trafo[trafo != 0]^(-1)
  }
  
  return(list(mat = mat, sig = sig))
}




#' Unique (rotated) IRF
#' 
#' Fixes permutations and signs of IRF as described in the article.
#'
#' @param theta doubles; deep parameters
#' @param tmpl SVARMA-WHF template, see \link{tmpl_whf_rev}
#' @param lag.max integer; maximal number of lags to be contained in the output
#'   Default set to 30.
#'
#' @return \link{pseries} object. IRF with 30 lags.
#' 
#' @export
irf_unique = function(theta, tmpl, lag.max = 30){
  mod = armamod_whf(theta, tmpl)
  ma = mod$polm_ma %r% mod$B
  
  ma0 = unclass(mod$polm_ma)[,,1]
  ma_normalized = ma %r% solve(ma0)
  
  out = unique_p_s(ma0)
  ma0new = out$mat %*% diag(out$sig)
  
  ma_new = ma_normalized %r% ma0new
  
  return(pseries(lmfd(a = mod$polm_ar, b = ma_new), lag.max = lag.max))
}


#' Permute or change sign of shocks to IRF
#'
#' First the permutation is right-multiplied on the IRF, then the sign matrix.
#' 
#' @param pseries pseries or array object
#' @param perm vector of dimension \code{dim_in}, describing the permutation to be postmultiplied on the transfer function
#' @param sign vector of dimension \code{dim_in} containing -1 and 1 as elements
#'
#' @return pseries object
#' @export
#'
#' @examples
#' ps = test_polm(dim = c(2,2), degree = 10, random = TRUE) %>% pseries()
#' irf_change(ps, c(2,1), c(1,-1))
irf_change = function(pseries, perm = NULL, sign = NULL){
  
  a = unclass(pseries)
  dim_a = dim(a)
  
  if (is.null(perm)){
    perm = 1:dim_a[2]
  }
  
  if (is.null(sign)){
    sign = rep(1L, dim_a[2])
  }
  
  stopifnot("*sign* must be a vector of dimension equal to the input dimension." = length(sign) == dim_a[2],
            "*perm* must be a vector of dimension equal to the input dimension." = length(perm) == dim_a[2],
            "Elements of *sign* must be equal to 1 or -1" = unique(sign) %in% c(-1,1),
            "Elements of *perm* need to be unique." = length(perm) == length(unique(perm)))
  
  a_list = purrr::map(1:dim_a[3], ~a[, perm, .x] %*% diag(sign) )
  
  for (ix in 1:dim_a[3]){
    a[,,ix] = a_list[[ix]]
  }
  
  a = structure(a, class = c('pseries','ratm'))
  return(a)
}

#' Skewness and kurtosis for SGT
#' 
#' Calculates skewness and kurtosis for the SGT family of distributions with mean 0 and variance 1 and parameters \code{l}, \code{p}, \code{q}.
#' Existence of moments depends on parameter values of \code{p}, \code{q}.
#' See \url{https://www.jstor.org/stable/2634700} and \url{https://cran.r-project.org/web/packages/sgt/vignettes/sgt.pdf}.
#'
#' @param l double between -1 and 1. Parametrises skewness
#' @param p,q positive double.  parametrise kurtosis
#'
#' @return List with slots \code{skewness} ,\code{kurtosis}.
#' @seealso \link[sgt]{sgt}
#' @export
#' 
#' @examples 
#' sgt_skew_kurt(-0.52,1.49,85616)
#' sgt_skew_kurt(0.15,1.92,7.89)
sgt_skew_kurt = function(l, p, q){
  
  # Factor for adjusting variance to 1
  v = q^(-1/p) *
    (
      (3*l^2+1) *
        (beta(3/p,q-2/p)/beta(1/p,q)) -
        4*l^2 *
        (beta(2/p,q-1/p)/beta(1/p,q))^2
    )^(-1/2)
  
  # See sgt vignette for the formulas
  skewness = 2 * q^(3/p) * l * v^3 / beta(1/p, q)^3 *
    (8*l^2*beta(2/p,q-1/p)^3 - 
       3*(1+3*l^2) * beta(1/p,q) * beta(2/p, q-1/p) * beta(3/p,q-2/p) + 
       2 * (1+l^2) * beta(1/p,q)^2 * beta(4/p, q-3/p)
    )
  
  kurtosis = q^(4/p) * v^4 / beta(1/p, q)^4 *
    (48*l^4*beta(2/p,q-1/p)^4 +
       24 * l^2 * (1+3*l^2) * beta(1/p,q) * beta(2/p, q-1/p)^2 * beta(3/p,q-2/p) -
       32 * l^2 * (1+l^2) * beta(1/p,q)^2 * beta(2/p, q-1/p) * beta(4/p,q-3/p) + 
       (1+10*l^2+5*l^4) * beta(1/p,q)^3 * beta(5/p, q-4/p)
    )
  
  return(list(skewness = skewness,
              kurtosis = kurtosis))
}
