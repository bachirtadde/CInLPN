#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <stdexcept>
#include <math.h>
#include <vector>
#include "genericfun.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// function that computes autocceralation between dimension from
// estimated parameters
//================================================================================
arma::mat autocorri(int K, arma::vec tau, arma::mat z0i, arma::mat zi, arma::mat matDw, 
arma::mat matDw_u, arma::mat matDu, arma::mat matB, arma::mat G_mat_A_0_to_tau_i,
arma::mat G_mat_prod_A_0_to_tau, double DeltaT){
  int m = tau.size();
  // found the max of the vector tau_i
  double maxTau = 0.0;
  for (unsigned i = 0; i < tau.size(); i++){
    if (tau[i] > maxTau)
      maxTau = tau[i];
  }
  int q = zi.n_cols;
  mat GrdZi = zeros((maxTau+1)*K, q);
  int p_j=0; // loop variable
  int p_k =0; // loop variable
  
  // ##### computering of GrdZi ####################################
  for(int t = 0; t<= maxTau; t++){
    if(t==0){
      GrdZi(span(t*K,(t+1)*K-1), span(0,q-1)) = 0*zi(span(t*K,(t+1)*K-1), span(0,q-1));
    }
    else{
      GrdZi(span(t*K,(t+1)*K-1), span(0,q-1)) = DeltaT*(zi(span(t*K,(t+1)*K-1), span(0,q-1)) +
      G_mat_A_0_to_tau_i(span(0,K-1),span(K*(t-1),t*K-1))*GrdZi(span((t-1)*K,t*K-1), span(0,q-1)));
    }
  }  
  // Computing of matVX_i=================================================================================
  mat matVX_i = zeros(K*m,K*m); // initialisation de matVX_i ? 0
  for( int j =0 ; j < m; j++){ 
    p_k = p_j;    
    for( int k =j ; k < m; k++){
      // ###### Computering of variance-covariance matrix the componante related to the RE ######## 
      vec vect = ones(K);
      mat prodA_0_to_t_j_1 = diagmat(vect);
      mat prodA_0_to_t_k_1 = diagmat(vect);
      if(tau(j)>0){
        prodA_0_to_t_j_1 = G_mat_prod_A_0_to_tau(span(0,K-1),span(K*(tau(j)-1),K*(tau(j)-1)+K-1));
      }
      if(tau(k)>0){
        prodA_0_to_t_k_1 = G_mat_prod_A_0_to_tau(span(0,K-1),span(K*(tau(k)-1),K*(tau(k)-1)+K-1));
      }
      matVX_i(span(p_j,(p_j+K-1)), span((p_k), (p_k+K-1))) +=
      pow(DeltaT, (tau(j)+tau(k)))*(prodA_0_to_t_j_1*z0i)*matDw*(z0i*prodA_0_to_t_k_1).t() + 
      pow(DeltaT, tau(j))*(prodA_0_to_t_j_1*z0i)*matDw_u*(GrdZi(span(tau(k)*K,(tau(k)+1)*K-1), span(0,q-1))).t() +
      pow(DeltaT, tau(k))*(GrdZi(span(tau(j)*K,(tau(j)+1)*K-1), span(0,q-1)))*matDw_u.t()*(z0i*prodA_0_to_t_k_1).t() +
      (GrdZi(span(tau(j)*K,(tau(j)+1)*K-1), span(0,q-1)))*matDu*(GrdZi(span(tau(k)*K,(tau(k)+1)*K-1), span(0,q-1))).t();
      //remplissage de la partie inf?reure de la matrice puisque matVX_i est sym?trique====
      if(p_j != p_k){
        matVX_i(span(p_k, (p_k+K-1)), span(p_j,(p_j+K-1))) = matVX_i(span(p_j,(p_j+K-1)), span(p_k, (p_k+K-1))).t();
      }
      p_k += K; // incr?mentation de p_k
    }
    p_j += K;// incr?mentation de p_j
  }
  //
  
  // calcul de l'autocorrelation================================== 
  //  // autocorrelation between latent processes
  mat corr_X = zeros(K*m,K*m);
  corr_X = diagmat(1/sqrt(matVX_i.diag()))*matVX_i*diagmat(1/sqrt(matVX_i.diag()));
  return(corr_X);
}



//==============================================================================================================================
/* for overall individuals

// [[Rcpp::export]]

*/

arma::mat autocorr(int K, arma::vec paraOpt, arma::vec paraFixe, arma::vec posfix, int ncol_x,
int ncol_x0, arma::mat zi, arma::vec q, int nb_paraD, arma::mat z0i, arma::vec q0,
arma::vec tau, arma::mat modA_mat_i, double DeltaT){
  
  int L = modA_mat_i.n_cols; // number of parameters per element of transition matrix A
  // Identifiacation of constraint parameters and non-constraint parameters
  int Nb_para = (int)posfix.size();
  vec paras = zeros(Nb_para);
  int Opt=0;
  int Fixe=0;
  for(int i=0; i < Nb_para ; i++){
    if(posfix[i]==0){
      paras[i] = paraOpt[Opt];
      Opt+=1;
    }
    else{
      paras[i] = paraFixe[Fixe];
      Fixe+=1;
    }
  }
  
  //Identification of groups of parameters
  int ipara =0;
  colvec alpha_mu0 = paras(span(ipara,ipara+ncol_x0-1));
  ipara += ncol_x0;
  colvec alpha_mu = DeltaT*paras(span(ipara,ipara+ncol_x-1));
  ipara += ncol_x;
  colvec alpha_D = paras(span(ipara,ipara + nb_paraD-1));
  ipara += nb_paraD;
  vec vec_alpha_ij = DeltaT*paras(span(ipara,ipara+L*K*K-1));
  ipara += L*K*K;
  vec paraB = zeros(K);
  vec paraSig = paras(span(ipara,ipara+K-1));
  ipara += K;
  //
  int nb_RE = sum(sum(q0)+sum(q));
  mat matD = DparChol(nb_RE, alpha_D);
  int n_cols_matD = matD.n_cols;
  mat matDw = matD(span(0,K-1),span(0,K-1));
  mat matDw_u = DeltaT*matD(span(0,K-1),span(K,n_cols_matD-1));
  mat matDu = DeltaT*matD(span(K,n_cols_matD-1),span(K,n_cols_matD-1))*DeltaT;
  mat matB = KmatDiag(paraB); // unstructured variance-covariance matrice
  mat Sig = KmatDiag(paraSig); // noice
  //Creation of matrix G_mat_prod_A_0_to_tau that contains all products  A(j) from t_i a Tmax: t_i \in 0, Tmax
  mat G_mat_prod_A_0_to_tau = GmatprodAstotau(K, vec_alpha_ij, tau, 0, DeltaT, modA_mat_i);
  mat G_mat_A_0_to_tau_i = GmatA0totaui(K, vec_alpha_ij, tau, DeltaT, modA_mat_i);
  //computering of individual corralations matrix 
  mat corr = autocorri(K, tau, z0i, zi, matDw, matDw_u, matDu, matB, G_mat_A_0_to_tau_i,
  G_mat_prod_A_0_to_tau, DeltaT);
  return(corr);
}
