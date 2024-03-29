// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// figtree
std::vector<double> figtree(std::vector<double> X, double h, std::vector<double> Q, std::vector<double> Y, double epsilon, std::vector<double> G);
RcppExport SEXP _TDApplied_figtree(SEXP XSEXP, SEXP hSEXP, SEXP QSEXP, SEXP YSEXP, SEXP epsilonSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Q(QSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(figtree(X, h, Q, Y, epsilon, G));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TDApplied_figtree", (DL_FUNC) &_TDApplied_figtree, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_TDApplied(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
