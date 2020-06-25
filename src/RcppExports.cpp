// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Loglik
double Loglik(int K, int nD, arma::vec& mapping, arma::vec& paraOpt, arma::vec& paraFixe, arma::vec& posfix, arma::vec& m_is, arma::mat& Mod_MatrixY, arma::mat& Mod_MatrixYprim, arma::vec& df, arma::mat& x, arma::mat& z, arma::vec& q, int nb_paraD, arma::mat& x0, arma::mat& z0, arma::vec& q0, arma::vec if_link, arma::vec& tau, arma::vec& tau_is, arma::mat& modA_mat, double DeltaT);
RcppExport SEXP _CInLPN_Loglik(SEXP KSEXP, SEXP nDSEXP, SEXP mappingSEXP, SEXP paraOptSEXP, SEXP paraFixeSEXP, SEXP posfixSEXP, SEXP m_isSEXP, SEXP Mod_MatrixYSEXP, SEXP Mod_MatrixYprimSEXP, SEXP dfSEXP, SEXP xSEXP, SEXP zSEXP, SEXP qSEXP, SEXP nb_paraDSEXP, SEXP x0SEXP, SEXP z0SEXP, SEXP q0SEXP, SEXP if_linkSEXP, SEXP tauSEXP, SEXP tau_isSEXP, SEXP modA_matSEXP, SEXP DeltaTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type nD(nDSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mapping(mappingSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type paraOpt(paraOptSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type paraFixe(paraFixeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type posfix(posfixSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m_is(m_isSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Mod_MatrixY(Mod_MatrixYSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Mod_MatrixYprim(Mod_MatrixYprimSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type nb_paraD(nb_paraDSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z0(z0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type q0(q0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type if_link(if_linkSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau_is(tau_isSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type modA_mat(modA_matSEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    rcpp_result_gen = Rcpp::wrap(Loglik(K, nD, mapping, paraOpt, paraFixe, posfix, m_is, Mod_MatrixY, Mod_MatrixYprim, df, x, z, q, nb_paraD, x0, z0, q0, if_link, tau, tau_is, modA_mat, DeltaT));
    return rcpp_result_gen;
END_RCPP
}
// fit
arma::mat fit(int K, int nD, arma::vec& mapping, arma::vec& paras, arma::vec& m_is, arma::mat& Mod_MatrixY, arma::vec df, arma::mat& x, arma::mat& z, arma::vec& q, int nb_paraD, arma::mat& x0, arma::mat& z0, arma::vec& q0, arma::vec if_link, arma::vec tau, arma::vec& tau_is, arma::mat& modA_mat, double DeltaT, int MCnr, arma::vec minY, arma::vec maxY, List& knots, arma::vec degree, double epsPred);
RcppExport SEXP _CInLPN_fit(SEXP KSEXP, SEXP nDSEXP, SEXP mappingSEXP, SEXP parasSEXP, SEXP m_isSEXP, SEXP Mod_MatrixYSEXP, SEXP dfSEXP, SEXP xSEXP, SEXP zSEXP, SEXP qSEXP, SEXP nb_paraDSEXP, SEXP x0SEXP, SEXP z0SEXP, SEXP q0SEXP, SEXP if_linkSEXP, SEXP tauSEXP, SEXP tau_isSEXP, SEXP modA_matSEXP, SEXP DeltaTSEXP, SEXP MCnrSEXP, SEXP minYSEXP, SEXP maxYSEXP, SEXP knotsSEXP, SEXP degreeSEXP, SEXP epsPredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type nD(nDSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mapping(mappingSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type paras(parasSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m_is(m_isSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Mod_MatrixY(Mod_MatrixYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type nb_paraD(nb_paraDSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z0(z0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type q0(q0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type if_link(if_linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau_is(tau_isSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type modA_mat(modA_matSEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    Rcpp::traits::input_parameter< int >::type MCnr(MCnrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type minY(minYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type maxY(maxYSEXP);
    Rcpp::traits::input_parameter< List& >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< double >::type epsPred(epsPredSEXP);
    rcpp_result_gen = Rcpp::wrap(fit(K, nD, mapping, paras, m_is, Mod_MatrixY, df, x, z, q, nb_paraD, x0, z0, q0, if_link, tau, tau_is, modA_mat, DeltaT, MCnr, minY, maxY, knots, degree, epsPred));
    return rcpp_result_gen;
END_RCPP
}
// pred
arma::mat pred(int K, int nD, arma::vec& mapping, arma::vec& paras, arma::vec& m_is, arma::mat& Mod_MatrixY, arma::vec df, arma::mat& x, arma::mat& z, arma::vec& q, int nb_paraD, arma::mat& x0, arma::mat& z0, arma::vec& q0, arma::vec if_link, arma::vec tau, arma::vec& tau_is, arma::mat& modA_mat, double DeltaT, int MCnr, arma::vec minY, arma::vec maxY, List& knots, arma::vec degree, double epsPred);
RcppExport SEXP _CInLPN_pred(SEXP KSEXP, SEXP nDSEXP, SEXP mappingSEXP, SEXP parasSEXP, SEXP m_isSEXP, SEXP Mod_MatrixYSEXP, SEXP dfSEXP, SEXP xSEXP, SEXP zSEXP, SEXP qSEXP, SEXP nb_paraDSEXP, SEXP x0SEXP, SEXP z0SEXP, SEXP q0SEXP, SEXP if_linkSEXP, SEXP tauSEXP, SEXP tau_isSEXP, SEXP modA_matSEXP, SEXP DeltaTSEXP, SEXP MCnrSEXP, SEXP minYSEXP, SEXP maxYSEXP, SEXP knotsSEXP, SEXP degreeSEXP, SEXP epsPredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type nD(nDSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mapping(mappingSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type paras(parasSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m_is(m_isSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Mod_MatrixY(Mod_MatrixYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type nb_paraD(nb_paraDSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z0(z0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type q0(q0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type if_link(if_linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau_is(tau_isSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type modA_mat(modA_matSEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    Rcpp::traits::input_parameter< int >::type MCnr(MCnrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type minY(minYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type maxY(maxYSEXP);
    Rcpp::traits::input_parameter< List& >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< double >::type epsPred(epsPredSEXP);
    rcpp_result_gen = Rcpp::wrap(pred(K, nD, mapping, paras, m_is, Mod_MatrixY, df, x, z, q, nb_paraD, x0, z0, q0, if_link, tau, tau_is, modA_mat, DeltaT, MCnr, minY, maxY, knots, degree, epsPred));
    return rcpp_result_gen;
END_RCPP
}
// vectorise
arma::vec vectorise(arma::mat& M);
RcppExport SEXP _CInLPN_vectorise(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(vectorise(M));
    return rcpp_result_gen;
END_RCPP
}
// KmatDiag
arma::mat KmatDiag(arma::vec& Kvector);
RcppExport SEXP _CInLPN_KmatDiag(SEXP KvectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type Kvector(KvectorSEXP);
    rcpp_result_gen = Rcpp::wrap(KmatDiag(Kvector));
    return rcpp_result_gen;
END_RCPP
}
// DparChol
arma::mat DparChol(int q, arma::vec& qvector);
RcppExport SEXP _CInLPN_DparChol(SEXP qSEXP, SEXP qvectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type qvector(qvectorSEXP);
    rcpp_result_gen = Rcpp::wrap(DparChol(q, qvector));
    return rcpp_result_gen;
END_RCPP
}
// vecaijt
arma::vec vecaijt(int K, int t, arma::vec& vec_alpha_ij, arma::mat& modA_mat);
RcppExport SEXP _CInLPN_vecaijt(SEXP KSEXP, SEXP tSEXP, SEXP vec_alpha_ijSEXP, SEXP modA_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vec_alpha_ij(vec_alpha_ijSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type modA_mat(modA_matSEXP);
    rcpp_result_gen = Rcpp::wrap(vecaijt(K, t, vec_alpha_ij, modA_mat));
    return rcpp_result_gen;
END_RCPP
}
// ConstrA
arma::mat ConstrA(int K, int t, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat);
RcppExport SEXP _CInLPN_ConstrA(SEXP KSEXP, SEXP tSEXP, SEXP DeltaTSEXP, SEXP vec_alpha_ijSEXP, SEXP modA_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vec_alpha_ij(vec_alpha_ijSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type modA_mat(modA_matSEXP);
    rcpp_result_gen = Rcpp::wrap(ConstrA(K, t, DeltaT, vec_alpha_ij, modA_mat));
    return rcpp_result_gen;
END_RCPP
}
// GmatA0totaui
arma::mat GmatA0totaui(int K, arma::vec& vec_alpha_ij, arma::vec& tau_i, double DeltaT, arma::mat modA_mat);
RcppExport SEXP _CInLPN_GmatA0totaui(SEXP KSEXP, SEXP vec_alpha_ijSEXP, SEXP tau_iSEXP, SEXP DeltaTSEXP, SEXP modA_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vec_alpha_ij(vec_alpha_ijSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau_i(tau_iSEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type modA_mat(modA_matSEXP);
    rcpp_result_gen = Rcpp::wrap(GmatA0totaui(K, vec_alpha_ij, tau_i, DeltaT, modA_mat));
    return rcpp_result_gen;
END_RCPP
}
// ProdA
arma::mat ProdA(int K, int t2, int t1, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat);
RcppExport SEXP _CInLPN_ProdA(SEXP KSEXP, SEXP t2SEXP, SEXP t1SEXP, SEXP DeltaTSEXP, SEXP vec_alpha_ijSEXP, SEXP modA_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< int >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vec_alpha_ij(vec_alpha_ijSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type modA_mat(modA_matSEXP);
    rcpp_result_gen = Rcpp::wrap(ProdA(K, t2, t1, DeltaT, vec_alpha_ij, modA_mat));
    return rcpp_result_gen;
END_RCPP
}
// GmatprodAstotau
arma::mat GmatprodAstotau(int K, arma::vec& vec_alpha_ij, arma::vec& tau, int t_ini, double DeltaT, arma::mat modA_mat);
RcppExport SEXP _CInLPN_GmatprodAstotau(SEXP KSEXP, SEXP vec_alpha_ijSEXP, SEXP tauSEXP, SEXP t_iniSEXP, SEXP DeltaTSEXP, SEXP modA_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vec_alpha_ij(vec_alpha_ijSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type t_ini(t_iniSEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type modA_mat(modA_matSEXP);
    rcpp_result_gen = Rcpp::wrap(GmatprodAstotau(K, vec_alpha_ij, tau, t_ini, DeltaT, modA_mat));
    return rcpp_result_gen;
END_RCPP
}
// tsGmatprodA0totau
arma::mat tsGmatprodA0totau(int K, arma::vec& vec_alpha_ij, arma::vec& tau, double DeltaT, arma::mat modA_mat);
RcppExport SEXP _CInLPN_tsGmatprodA0totau(SEXP KSEXP, SEXP vec_alpha_ijSEXP, SEXP tauSEXP, SEXP DeltaTSEXP, SEXP modA_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vec_alpha_ij(vec_alpha_ijSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type modA_mat(modA_matSEXP);
    rcpp_result_gen = Rcpp::wrap(tsGmatprodA0totau(K, vec_alpha_ij, tau, DeltaT, modA_mat));
    return rcpp_result_gen;
END_RCPP
}
// compoYiNA
arma::vec compoYiNA(arma::mat& Yi);
RcppExport SEXP _CInLPN_compoYiNA(SEXP YiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Yi(YiSEXP);
    rcpp_result_gen = Rcpp::wrap(compoYiNA(Yi));
    return rcpp_result_gen;
END_RCPP
}
// YiwoNA
arma::vec YiwoNA(arma::vec Yi);
RcppExport SEXP _CInLPN_YiwoNA(SEXP YiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Yi(YiSEXP);
    rcpp_result_gen = Rcpp::wrap(YiwoNA(Yi));
    return rcpp_result_gen;
END_RCPP
}
// matNui
arma::mat matNui(int nD, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0, arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i);
RcppExport SEXP _CInLPN_matNui(SEXP nDSEXP, SEXP tau_iSEXP, SEXP DeltaTSEXP, SEXP x0iSEXP, SEXP alpha_mu0SEXP, SEXP xiSEXP, SEXP alpha_muSEXP, SEXP G_mat_A_0_to_tau_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nD(nDSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau_i(tau_iSEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x0i(x0iSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha_mu0(alpha_mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha_mu(alpha_muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G_mat_A_0_to_tau_i(G_mat_A_0_to_tau_iSEXP);
    rcpp_result_gen = Rcpp::wrap(matNui(nD, tau_i, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i));
    return rcpp_result_gen;
END_RCPP
}
// f_Yi_r_NA_by0
arma::mat f_Yi_r_NA_by0(arma::mat& Yi);
RcppExport SEXP _CInLPN_f_Yi_r_NA_by0(SEXP YiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Yi(YiSEXP);
    rcpp_result_gen = Rcpp::wrap(f_Yi_r_NA_by0(Yi));
    return rcpp_result_gen;
END_RCPP
}
// YiNui
arma::vec YiNui(int nD, arma::mat matrixP, arma::vec& tau, arma::vec& tau_i, double DeltaT, arma::mat& Yi, arma::mat& x0i, arma::colvec& alpha_mu0, arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i);
RcppExport SEXP _CInLPN_YiNui(SEXP nDSEXP, SEXP matrixPSEXP, SEXP tauSEXP, SEXP tau_iSEXP, SEXP DeltaTSEXP, SEXP YiSEXP, SEXP x0iSEXP, SEXP alpha_mu0SEXP, SEXP xiSEXP, SEXP alpha_muSEXP, SEXP G_mat_A_0_to_tau_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nD(nDSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type matrixP(matrixPSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau_i(tau_iSEXP);
    Rcpp::traits::input_parameter< double >::type DeltaT(DeltaTSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Yi(YiSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x0i(x0iSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha_mu0(alpha_mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha_mu(alpha_muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G_mat_A_0_to_tau_i(G_mat_A_0_to_tau_iSEXP);
    rcpp_result_gen = Rcpp::wrap(YiNui(nD, matrixP, tau, tau_i, DeltaT, Yi, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i));
    return rcpp_result_gen;
END_RCPP
}
// VecToMat
arma::mat VecToMat(arma::vec& y, int K, int m_i);
RcppExport SEXP _CInLPN_VecToMat(SEXP ySEXP, SEXP KSEXP, SEXP m_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type m_i(m_iSEXP);
    rcpp_result_gen = Rcpp::wrap(VecToMat(y, K, m_i));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CInLPN_Loglik", (DL_FUNC) &_CInLPN_Loglik, 22},
    {"_CInLPN_fit", (DL_FUNC) &_CInLPN_fit, 25},
    {"_CInLPN_pred", (DL_FUNC) &_CInLPN_pred, 25},
    {"_CInLPN_vectorise", (DL_FUNC) &_CInLPN_vectorise, 1},
    {"_CInLPN_KmatDiag", (DL_FUNC) &_CInLPN_KmatDiag, 1},
    {"_CInLPN_DparChol", (DL_FUNC) &_CInLPN_DparChol, 2},
    {"_CInLPN_vecaijt", (DL_FUNC) &_CInLPN_vecaijt, 4},
    {"_CInLPN_ConstrA", (DL_FUNC) &_CInLPN_ConstrA, 5},
    {"_CInLPN_GmatA0totaui", (DL_FUNC) &_CInLPN_GmatA0totaui, 5},
    {"_CInLPN_ProdA", (DL_FUNC) &_CInLPN_ProdA, 6},
    {"_CInLPN_GmatprodAstotau", (DL_FUNC) &_CInLPN_GmatprodAstotau, 6},
    {"_CInLPN_tsGmatprodA0totau", (DL_FUNC) &_CInLPN_tsGmatprodA0totau, 5},
    {"_CInLPN_compoYiNA", (DL_FUNC) &_CInLPN_compoYiNA, 1},
    {"_CInLPN_YiwoNA", (DL_FUNC) &_CInLPN_YiwoNA, 1},
    {"_CInLPN_matNui", (DL_FUNC) &_CInLPN_matNui, 8},
    {"_CInLPN_f_Yi_r_NA_by0", (DL_FUNC) &_CInLPN_f_Yi_r_NA_by0, 1},
    {"_CInLPN_YiNui", (DL_FUNC) &_CInLPN_YiNui, 11},
    {"_CInLPN_VecToMat", (DL_FUNC) &_CInLPN_VecToMat, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_CInLPN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
