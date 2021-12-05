// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_svarmawhf_RCPPEXPORTS_H_GEN_
#define RCPP_svarmawhf_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace svarmawhf {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("svarmawhf", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("svarmawhf", "_svarmawhf_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in svarmawhf");
            }
        }
    }

    inline arma::mat get_residuals(const arma::mat& data_in, const arma::mat& polm_ar, arma::mat& polm_ma_bwd, const arma::mat& polm_ma_fwd, const arma::uword& kappa, const arma::uword& k) {
        typedef SEXP(*Ptr_get_residuals)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_get_residuals p_get_residuals = NULL;
        if (p_get_residuals == NULL) {
            validateSignature("arma::mat(*get_residuals)(const arma::mat&,const arma::mat&,arma::mat&,const arma::mat&,const arma::uword&,const arma::uword&)");
            p_get_residuals = (Ptr_get_residuals)R_GetCCallable("svarmawhf", "_svarmawhf_get_residuals");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_get_residuals(Shield<SEXP>(Rcpp::wrap(data_in)), Shield<SEXP>(Rcpp::wrap(polm_ar)), Shield<SEXP>(Rcpp::wrap(polm_ma_bwd)), Shield<SEXP>(Rcpp::wrap(polm_ma_fwd)), Shield<SEXP>(Rcpp::wrap(kappa)), Shield<SEXP>(Rcpp::wrap(k)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline void solve_ARMA_cpp(const arma::mat& a, const arma::mat& b, const arma::mat& u, arma::mat& y, int t0) {
        typedef SEXP(*Ptr_solve_ARMA_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_solve_ARMA_cpp p_solve_ARMA_cpp = NULL;
        if (p_solve_ARMA_cpp == NULL) {
            validateSignature("void(*solve_ARMA_cpp)(const arma::mat&,const arma::mat&,const arma::mat&,arma::mat&,int)");
            p_solve_ARMA_cpp = (Ptr_solve_ARMA_cpp)R_GetCCallable("svarmawhf", "_svarmawhf_solve_ARMA_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_solve_ARMA_cpp(Shield<SEXP>(Rcpp::wrap(a)), Shield<SEXP>(Rcpp::wrap(b)), Shield<SEXP>(Rcpp::wrap(u)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(t0)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
    }

}

#endif // RCPP_svarmawhf_RCPPEXPORTS_H_GEN_