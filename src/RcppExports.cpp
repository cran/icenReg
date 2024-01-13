// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ic_parList
Rcpp::List ic_parList(Rcpp::List R_parList);
RcppExport SEXP _icenReg_ic_parList(SEXP R_parListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type R_parList(R_parListSEXP);
    rcpp_result_gen = Rcpp::wrap(ic_parList(R_parList));
    return rcpp_result_gen;
END_RCPP
}
// R_ic_bayes
Rcpp::List R_ic_bayes(Rcpp::List R_bayesList, Rcpp::Function priorFxn, Rcpp::List R_ic_parList);
RcppExport SEXP _icenReg_R_ic_bayes(SEXP R_bayesListSEXP, SEXP priorFxnSEXP, SEXP R_ic_parListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type R_bayesList(R_bayesListSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type priorFxn(priorFxnSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type R_ic_parList(R_ic_parListSEXP);
    rcpp_result_gen = Rcpp::wrap(R_ic_bayes(R_bayesList, priorFxn, R_ic_parList));
    return rcpp_result_gen;
END_RCPP
}
// computeConditional_p
Rcpp::NumericVector computeConditional_p(Rcpp::NumericVector q, Rcpp::NumericVector etas, Rcpp::NumericMatrix baselineParams, Rcpp::CharacterVector reg_model, Rcpp::CharacterVector base_dist);
RcppExport SEXP _icenReg_computeConditional_p(SEXP qSEXP, SEXP etasSEXP, SEXP baselineParamsSEXP, SEXP reg_modelSEXP, SEXP base_distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type etas(etasSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type baselineParams(baselineParamsSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type reg_model(reg_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type base_dist(base_distSEXP);
    rcpp_result_gen = Rcpp::wrap(computeConditional_p(q, etas, baselineParams, reg_model, base_dist));
    return rcpp_result_gen;
END_RCPP
}
// computeConditional_q
Rcpp::NumericVector computeConditional_q(Rcpp::NumericVector p, Rcpp::NumericVector etas, Rcpp::NumericMatrix baselineParams, Rcpp::CharacterVector reg_model, Rcpp::CharacterVector base_dist);
RcppExport SEXP _icenReg_computeConditional_q(SEXP pSEXP, SEXP etasSEXP, SEXP baselineParamsSEXP, SEXP reg_modelSEXP, SEXP base_distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type etas(etasSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type baselineParams(baselineParamsSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type reg_model(reg_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type base_dist(base_distSEXP);
    rcpp_result_gen = Rcpp::wrap(computeConditional_q(p, etas, baselineParams, reg_model, base_dist));
    return rcpp_result_gen;
END_RCPP
}
