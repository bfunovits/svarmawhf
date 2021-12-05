// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]

// To make your C++ code callable from C++ code in other packages.
// This will generate a header file, inst/include/mypackage.h that
// can be included by other packages
// [[Rcpp::interfaces(r, cpp)]]


// via the exports attribute we tell Rcpp to make this function
// available from R

//' Rcpp version of solve_ARMA
//'
//' @param a,b Parameter matrices. Check R version for details.
//' @param u,y Data. y is overwritten, u is const reference. 
//' @param t0 Integer specifying at what index iteration should be started.
//' 
//' @return void
//'
//' @export
//' @name solve_ARMA
//' @keywords internal
// [[Rcpp::export]]
void solve_ARMA_cpp(const arma::mat& a, const arma::mat& b, 
                    const arma::mat& u, arma::mat& y, 
                    int t0) {
  int m = y.n_rows;
  int n = u.n_rows;
  int nobs = y.n_cols;

  int p = a.n_cols / m;
  int q = -1;
  if (n > 0) {
    q = (b.n_cols / n) - 1;
  }
  // Rcout << "m: " << m << " n: " << n << " p: " << p << " q: " << q << " nobs: " << nobs << std::endl;

  int shift = 0;

  // vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
  arma::vec yvec = arma::vec(y.memptr(), m*nobs, false, true);

  // y[t] = b[0] u[t]
  // Rcout << "t0-1: " << t0-1 << std::endl;
  if ( (n*(q+1)) > 0 ) {
    y.cols(t0-1, nobs-1) = b.cols(0, n-1) * u.cols(t0-1, nobs-1);
    for (int i = 1; i <= q; i++) {
      // y[t] = y[t] + b[i] u[t-i]
      // take care of missing initial values
      shift = std::max(i+1-t0,0);
      // Rcout << "shift+t0-1: " << shift+t0-1 << std::endl;
      y.cols(shift+t0-1, nobs-1) = y.cols(shift+t0-1, nobs-1) + b.cols(i*n, (i+1)*n-1) * u.cols(shift+t0-1-i, nobs-1-i);
    }
  } else {
    // Rcout << "t0-1: " << t0-1 << " nobs-1: " << nobs-1 << " ncols(y): " << y.n_cols << " nrows(y): " << y.n_rows << std::endl;
    y.cols(t0-1, nobs-1).zeros();
  }

  if (p > 0) {
    int j1 = (t0-1-p)*m;
    int j2 = (t0-1)*m-1;
    int i1 = (t0-1)*m;
    int i2 = t0*m-1;
    int t = t0 - 1;
    while (t < p) {
      // y[t] = u[t] + (a[p],...,a[1])(y[t-p]',...,y[t-1]')'
      // take care of missing initial values
      shift = std::max(-j1,0);
      // Rcout << "t: " << t << " shift: " << shift << " j1:j2"  << j1 << ":"  << j2 << " i1:i2 "  << i1 << ":"  << i2 << std::endl;
      if (j2 >= 0) {
        yvec.subvec(i1, i2) = yvec.subvec(i1, i2) + a.cols(shift, m*p-1) * yvec.subvec(0, j2);
      }

      t++;
      j1 = j1 + m;
      j2 = j2 + m;
      i1 = i1 + m;
      i2 = i2 + m;
    }
    while (t < nobs) {
      // y[t] = u[t] + (a[p],...,a[1])(y[t-p]',...,y[t-1]')'
      // Rcout << "t: " << t << " "  << j1 << " "  << j2 << " "  << i1 << " "  << i2 << std::endl;
      yvec.subvec(i1, i2) = yvec.subvec(i1, i2) + a * yvec.subvec(j1, j2);

      t++;
      j1 = j1 + m;
      j2 = j2 + m;
      i1 = i1 + m;
      i2 = i2 + m;
    }
  }
  return;
}

