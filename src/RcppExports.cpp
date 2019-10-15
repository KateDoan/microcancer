// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/microcancer.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// create_IMABC
void create_IMABC();
static SEXP _microcancer_create_IMABC_try() {
BEGIN_RCPP
    create_IMABC();
    return R_NilValue;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _microcancer_create_IMABC() {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_microcancer_create_IMABC_try());
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// check_sim
std::vector<double> check_sim();
static SEXP _microcancer_check_sim_try() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    rcpp_result_gen = Rcpp::wrap(check_sim());
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _microcancer_check_sim() {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_microcancer_check_sim_try());
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// squared_dist_vec
double squared_dist_vec(std::vector<double> v1, std::vector<double> v2);
RcppExport SEXP _microcancer_squared_dist_vec(SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(squared_dist_vec(v1, v2));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm_arma
arma::vec dmvnrm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd);
RcppExport SEXP _microcancer_dmvnrm_arma(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_arma(x, mean, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _microcancer_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("void(*create_IMABC)()");
        signatures.insert("std::vector<double>(*check_sim)()");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _microcancer_RcppExport_registerCCallable() { 
    R_RegisterCCallable("microcancer", "_microcancer_create_IMABC", (DL_FUNC)_microcancer_create_IMABC_try);
    R_RegisterCCallable("microcancer", "_microcancer_check_sim", (DL_FUNC)_microcancer_check_sim_try);
    R_RegisterCCallable("microcancer", "_microcancer_RcppExport_validate", (DL_FUNC)_microcancer_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_microcancer_create_IMABC", (DL_FUNC) &_microcancer_create_IMABC, 0},
    {"_microcancer_check_sim", (DL_FUNC) &_microcancer_check_sim, 0},
    {"_microcancer_squared_dist_vec", (DL_FUNC) &_microcancer_squared_dist_vec, 2},
    {"_microcancer_dmvnrm_arma", (DL_FUNC) &_microcancer_dmvnrm_arma, 4},
    {"_microcancer_RcppExport_registerCCallable", (DL_FUNC) &_microcancer_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_microcancer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
