#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]


//' Obtain residuals with Rcpp
//' 
//' For given AR, bwd MA, fwd MA parameters, calculate the residuals (B and Sigma are taken into account later).
//' The input arguments \code{polm_ar}, \code{polm_ma_fwd}, \code{polm_ma_bwd} are in a different format than usual for this Rcpp function, 
//' i.e. they are in wide matrices.
//' Also note that \eqn{p_0} is not necessarily the identity matrix and needs to be inverted in the code below.
//' \cr
//' What is happening in this function is better documented in its R version: \link{get_residuals_R}
//'
//' @param data_in const arma::mat reference. Observations correspond to columns.
//' @param polm_ar,polm_ma_fwd const arma::mat references. Wide matrices.
//'   \eqn{I_n, a_1, \ldots, a_p} and \eqn{f_1, \ldots, f_K}
//' @param polm_ma_bwd arma::mat reference. 
//'   \strong{NOT CONST} because \eqn{p_0} needs to be inverted, 
//'   applied to \eqn{p_q, \ldots, p_1, p_0}, and then the new \eqn{p_0 = I_n} needs to be discarded
//' @param kappa,k const arma::uword references. Partial indices
//' 
//' @return Object of type arma::mat. Wide data matrix (rows correspond to variables, observations to columns) of residuals.
//'
//' @export
// [[Rcpp::export]]
arma::mat get_residuals(const arma::mat& data_in, 
                        const arma::mat& polm_ar, arma::mat& polm_ma_bwd, const arma::mat& polm_ma_fwd, 
                        const arma::uword& kappa, const arma::uword& k) {
  
  arma::uword dim_out = data_in.n_rows;
  arma::uword n_obs = data_in.n_cols;
  
  arma::mat data_out = arma::zeros(dim_out, n_obs);
  
      // Rcpp::Rcerr << "dim_out: " << dim_out << " n_obs: " << n_obs << std::endl;
  
  arma::uword ARorder = polm_ar.n_cols/dim_out - 1;
  arma::uword MAorder_bwd = polm_ma_bwd.n_cols/dim_out - 1;
  arma::uword MAorder_fwd = polm_ma_fwd.n_cols/dim_out;
  
  // AR coefficient in wide matrix form: (I, a_1, ..., a_p)
  // First *ARorder* observations in data_out remain zero -> should not be used anymore
  for (arma::uword ix_coeff = 0; ix_coeff <= ARorder; ix_coeff++) {
    data_out.cols(ARorder, n_obs-1) = data_out.cols(ARorder, n_obs-1) + 
      polm_ar.cols(ix_coeff*dim_out, (ix_coeff+1)*dim_out-1) * data_in.cols(ARorder-ix_coeff, (n_obs-1-ix_coeff));
  }
  
  // Rcpp::Rcerr << "Everything ok after first AR loop" << std::endl;
  // Rcpp::Rcerr << data_out << std::endl;
  
  
  // p_0 not identity!!!!
  // using solve() and trimatl(), to indicate that the matrix is lower-triangular, is in general faster and more accurate than using inv() or the member function .i()
  // see http://arma.sourceforge.net/docs.html#i_member and http://arma.sourceforge.net/docs.html#solve
  // However, I am unsure if immediate reassignment is fine here. Thus I still use inv()
  
  // Obtain p_0 and premultiply inv(p_0) on data and (p_q, p_{q-1}, ... , p_1)
  arma::mat p0inv = polm_ma_bwd.tail_cols(dim_out).i();
  
  // Rcpp::Rcerr << p0inv << std::endl;
  
  data_out = p0inv * data_out;
  polm_ma_bwd = p0inv * polm_ma_bwd;
  
  // data_out = solve(p0, data_out);
  // polm_ma_bwd = solve(p0, polm_ma_bwd);
    
  // Discard new p_0 (which is equal to the identity matrix) from new backward MA poly (p_q, p_{q-1}, ... , p_1, I)
  arma::mat polm_ma_bwd_no0 = polm_ma_bwd.head_cols(dim_out * MAorder_bwd);
  
  //Rcpp::Rcerr << polm_ma_bwd_no0 << std::endl;
  
  
  // Vectorize (wide) data matrix with shared memory
  // Using arma::vectorise() in the loops below would have been easier...
  
  // arma::vec vec_data_out = arma::vec(data_out.memptr(), dim_out * n_obs, false, true);
  
  // MA bwd coefficient need to be reversed, no zero-lag coefficient matrix: (p_q, p_{q-1}, ... , p_1)
  if (MAorder_bwd>0){
    
    for (arma::uword ix_t = ARorder+1; ix_t < ARorder+MAorder_bwd; ix_t++){
      data_out.col(ix_t) =
        data_out.col(ix_t) - 
        polm_ma_bwd_no0.tail_cols(dim_out*(ix_t-ARorder)) * vectorise(data_out.cols(ARorder, ix_t-1));
        // vec_data_out.subvec(dim_out*ARorder, dim_out*ix_t-1);
    }

        // Rcpp::Rcerr << "Everything ok after first MAbwd loop" << std::endl;

    for (arma::uword ix_t = ARorder+MAorder_bwd; ix_t < n_obs; ix_t++){
      data_out.col(ix_t) = data_out.col(ix_t) - polm_ma_bwd_no0 * vectorise(data_out.cols(ix_t-MAorder_bwd, ix_t-1));
      }
    
    // Rcpp::Rcerr << "Everything ok after second MAbwd loop" << std::endl;
    // Rcpp::Rcerr << data_out << std::endl;
    
  }
  
  // Apply shift if k > 0 (kappa has no effect on programming/Toeplitz calculations)
  if (k>0){
    arma::mat shifted = arma::shift(data_out.rows(k,dim_out-1), 1, 1);
    data_out.rows(k,dim_out-1) = shifted;
    data_out.col(0).zeros();
  }

  // MA fwd coefficients DO NOT need to be reversed, but also no zero-lag coefficient matrix: (f_1, f_2, ... , f_q)
  if (MAorder_fwd>0){
    
    for (arma::uword ix_t = n_obs-1-1; ix_t > n_obs-1-MAorder_fwd; ix_t--){
      
      // Rcpp::Rcerr << "ix_t: " << ix_t << std::endl;
      // Rcpp::Rcerr << "Poly: " << polm_ma_fwd.head_cols(dim_out*(n_obs-1-ix_t)) << std::endl;
      // Rcpp::Rcerr << "Vec: " << vec_data_out.subvec(dim_out*(ix_t+1), dim_out*(n_obs)-1) << std::endl;
      // Rcpp::Rcerr << "Index start: " << dim_out*(ix_t+1) << std::endl;
      // Rcpp::Rcerr << "Index end: " << dim_out*(n_obs)-1 << std::endl;
      
      data_out.col(ix_t) = data_out.col(ix_t) - polm_ma_fwd.head_cols(dim_out*(n_obs-1-ix_t)) * vectorise(data_out.cols(ix_t+1, n_obs-1));
        // vec_data_out.subvec(dim_out*(ix_t+1), dim_out*(n_obs)-1);
    }

    // Rcpp::Rcerr << "Everything ok after first MAfwd loop" << std::endl;
    // Rcpp::Rcerr << data_out << std::endl;
    
    for (arma::uword ix_t = n_obs-1-MAorder_fwd; ix_t > 0; ix_t--){
      data_out.col(ix_t) = data_out.col(ix_t) - polm_ma_fwd * data_out.cols(ix_t+1, ix_t+MAorder_fwd).as_col();
        // vec_data_out.subvec( dim_out*(ix_t+1), dim_out*(ix_t+1+MAorder_fwd)-1);
    }
    
    // Rcpp::Rcerr << "Everything ok after second MAfwd loop" << std::endl;
    // Rcpp::Rcerr << data_out << std::endl;
    
  }
  
  return data_out;
  
}
