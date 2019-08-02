// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// simA1
SEXP simA1(const double Q, const Eigen::Map<Eigen::VectorXd> D, const Eigen::Map<Eigen::VectorXd> X, const Eigen::Map<Eigen::VectorXd> Tx, const Eigen::Map<Eigen::VectorXd> inputA);
RcppExport SEXP _pr2ar_simA1(SEXP QSEXP, SEXP DSEXP, SEXP XSEXP, SEXP TxSEXP, SEXP inputASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type D(DSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type Tx(TxSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type inputA(inputASEXP);
    rcpp_result_gen = Rcpp::wrap(simA1(Q, D, X, Tx, inputA));
    return rcpp_result_gen;
END_RCPP
}
// simA2
SEXP simA2(const double Q, const Eigen::Map<Eigen::VectorXd> D, const Eigen::Map<Eigen::VectorXd> X, const Eigen::Map<Eigen::VectorXd> Tx, const Eigen::Map<Eigen::VectorXd> inputA);
RcppExport SEXP _pr2ar_simA2(SEXP QSEXP, SEXP DSEXP, SEXP XSEXP, SEXP TxSEXP, SEXP inputASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type D(DSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type Tx(TxSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type inputA(inputASEXP);
    rcpp_result_gen = Rcpp::wrap(simA2(Q, D, X, Tx, inputA));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pr2ar_simA1", (DL_FUNC) &_pr2ar_simA1, 5},
    {"_pr2ar_simA2", (DL_FUNC) &_pr2ar_simA2, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_pr2ar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
