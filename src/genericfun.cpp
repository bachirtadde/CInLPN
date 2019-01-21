// #define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <stdexcept>
#include <math.h>
#include <stdlib.h> /* srand, rand */
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

//===========================================================================================================
//function print matrix
int f_mat_print( arma::mat& B)
{
  int J = B.n_cols;
  int I = B.n_rows;
  printf("\n \n===========================\n affichage de la matrice\n============================\n \n");
  printf(" dimension ::::%d * %d\n\n", I,J);
  for(int i =0; i<I; i++)
  {
    printf("| ");
    for(int j =0; j<J; j++)
    {
      printf("%.4f ",B(i,j));
    }
    printf("| \n");
  }
  printf("\n \n===========================\n fin de l'affichage\n============================\n \n");
  return 1;
}

//===========================================================================================================
// function to test if a matrix is inversible
arma::mat f_inv_mat(arma::mat& B)
{
  int g = 0;
  mat C;
  try{

    C = inv_sympd(B);
    //C = inv(B);

  }
  catch(const std::runtime_error& e){
    g=1;
    printf("Matrix B is not inversible : %d \n", g);
    C = zeros(B.n_cols,B.n_cols);
  }
  //   f_mat_print(C);
  if(g==1){
    printf("Matrix B is not inversible : %d \n", g);
  }

  return (C);
}


//======================================================================
//' Function that vectorises a matrix by rows
//' 
//' @param M a matrix
//' @export
// [[Rcpp::export]]
arma::vec vectorise(arma::mat& M){
  int n_c = M.n_cols;
  int n_r = M.n_rows;
  vec v = zeros(n_r*n_c);
  int pp=0;
  for(int i =0; i<n_r; i++){
    for(int k =0; k<n_c; k++){
      v(pp) = M(i,k);
      pp++;
    }
  }
  return(v);
}

//===========================================================================================================
arma::vec InnerProd(arma::vec v1, arma::vec v2){
  int n = v1.size();
  vec ip = zeros(n);
  for(int i=0; i < n; i++){
    ip[i] = v1[i]*v2[i];
  }
  return (ip);
}


//===========================================================================================================
//' Function that creates a K-block diagonal matrix with non diagonal element fixed to 0
//' 
//' @param Kvector a vector of length K
//' 
//' @return a diagonal matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat KmatDiag(arma::vec& Kvector){
  // Kvector : K-vector of paramters used to construct the matrix Kmat
  for(int i=0; i < (int)Kvector.size(); i++){
    Kvector[i] = pow(Kvector[i],2);
  }
  mat Kmat = diagmat(Kvector);
  return (Kmat);
}


//'  Function that computes a symetric D matric from it Cholesky L
//'  
//' @param q an integer
//' @param qvector a vector of length q
//' 
//' @return a symetric positive define matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat DparChol(int q, arma::vec& qvector){
  // qvector : vecteur de param?tres de taille K pour la construction de la matrice
  mat L = zeros<mat>(q,q);
  int pc=0 ;//  pour index? un param?tre dans le vecteur paraChol
  for( int j=0;j< q; ++j){
    for(int i=j; i<q; ++i){
      L(i,j) = qvector[pc];
      pc+=1;
    }
  }
  mat D = L*L.t();
  return (D);
}

//===========================================================================================================
/* ***********************************
Function aijt
Description:
aijt a function that computes the coefficient a_ij of the transition matrix A at time t from
model.matrix and vector of parameters
*/

double aijt(int t, arma::vec alpha_ijl,  arma::mat modA_mat){
  //modA_mat: model.matrix
  // remember that, length(alpha_ijl) = ncol(modA_mat)
  double a_ij_t = as_scalar(modA_mat(span(t,t),span(0,modA_mat.n_cols-1))*alpha_ijl);
  return(a_ij_t);
}



//'  Function that computes all coefficients of the transition matrix at time t
//'  
//' @param K an integer representing the size of K*K matrix
//' @param t an integer indicating the time at which coefficients are computed 
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return a vector of elements of the transition matrix
//' @export
//' 
// [[Rcpp::export]]
arma::vec vecaijt( int K, int t, arma::vec& vec_alpha_ij, arma::mat& modA_mat){
  // K = number of outcomes
  // L = taille de la matrice modA_mat.
  int L = modA_mat.n_cols;
  arma::vec vec_a_ij_t = zeros<vec>(K*K);
  int pp = 0; // index of loop
  for( int k = 0 ; k < K*K; k++){
    vec_a_ij_t(k) = aijt(t, vec_alpha_ij(span(pp,pp+L-1)), modA_mat);
    pp+=L;
  }
  return(vec_a_ij_t);

}

//===========================================================================================================

//'  Function that  constructs the transition matrix at time t
//'  
//' @param K an integer representing the size of K*K matrix
//' @param t an integer indicating the time at which coefficients are computed 
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return Id + DeltaT*A where A is the transition matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat ConstrA(int K, int t, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat){
  vec a_ij_t = vecaijt(K, t, vec_alpha_ij, modA_mat);
  mat A = zeros<mat>(K,K);
  vec vect = ones(K);
  mat Id = diagmat(vect);
  int p=0;
  for(int i = 0;i< K;i++){
    for(int j = 0;j< K;j++){
      A(i,j) = a_ij_t(p);
      p++;
    }
  }
  return(Id + DeltaT*A); // Atild = ((1/DeltaT)*Id + A)
}

//=======================================================================================================

//'  Function that  creates a matrix K,K*(max(tau_i)-1) containing  sub-matrices {(A_t)}_{t=0,tau_i-1}
//'  
//' @param K an integer representing the size of K*K matrix
//' @param tau_i vector of integers indicating times at which coefficients are computed 
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return a matrix containing matrix of form Id + DeltaT*A
//' @export
//' 
// [[Rcpp::export]]
arma::mat GmatA0totaui(int K, arma::vec& vec_alpha_ij, arma::vec& tau_i, double DeltaT,  arma::mat modA_mat) {
  int siz_tau_i = tau_i.size();

  mat G_mat_A_0_to_tau_i = zeros<mat>(K,K*siz_tau_i);
  for (int i=0; i< siz_tau_i; i++) {
    G_mat_A_0_to_tau_i(span(0,K-1),span(K*i,K*i+K-1)) = ConstrA(K, tau_i[i], DeltaT, vec_alpha_ij, modA_mat);
  }
  return (G_mat_A_0_to_tau_i);
}

//=======================================================================================================

//'  Function that  computes the product of A(t) for t1 to t2
//'  
//' @param K an integer representing the size of K*K matrix
//' @param t1 indicates the started discretized time
//' @param t2 indicates the ended discretized time
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return a matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat ProdA(int K, int t2, int t1, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat){
  // t2 = end time
  // t1 = start time
  // t1 <= t2
  mat pA = ConstrA(K,t1, DeltaT, vec_alpha_ij, modA_mat);
  for(int i = (t1+1); i <= t2; i++){
    pA = pA*ConstrA(K, i, DeltaT, vec_alpha_ij, modA_mat);
  }
  return(pA);
}


//===========================================================================================
//'  Function that  creates a big matrix containing  Prod(A_t)t=t_ini,tau.
//'  
//' @param K an integer representing the size of K*K matrix
//' @param t_ini indicates the started discretized time
//' @param tau vector of integers indicating times 
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return a matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat GmatprodAstotau( int K, arma::vec& vec_alpha_ij, arma::vec& tau,
                           int t_ini, double DeltaT, arma::mat modA_mat) {
  // t_ini = t initial
  int nb_t = tau.size()-t_ini;
  int ii = 0; //index de boucle
  mat G_mat_prod_A_s_to_tau = zeros<mat>(K,K*nb_t);
  for (int i=t_ini; i< (int)tau.size(); i++) {
    G_mat_prod_A_s_to_tau(span(0,(K-1)),span(ii,(ii+K-1))) = ProdA(K, tau[i], t_ini, DeltaT, vec_alpha_ij, modA_mat);
    ii +=K;
  }
  return (G_mat_prod_A_s_to_tau);
}

//===========================================================================================
//'  Function that  creates a big matrix ts_G_mat_prod_A_0_to_tau containing  Prod(A_t)t=0,tau.
//'  
//' @param K an integer representing the size of K*K matrix
//' @param tau vector of integers indicating times 
//' @param DeltaT double that indicates the discretization step
//' @param modA_mat model.matrix for elements of the transistion matrix  
//' @param vec_alpha_ij a vector of overall parameters associated to the
//' model.matrix for elements of the transistion matrix 
//' 
//' @return a matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat tsGmatprodA0totau(int K, arma::vec& vec_alpha_ij, arma::vec& tau, double DeltaT, arma::mat modA_mat) {
  //Creation of matrix G_mat_prod_A_0_to_tau that contains all products  A(j) from t_j a Tmax: t_j \in 0, Tmax
  int T = tau.size();
  mat ts_G_mat_prod_A_0_to_tau = zeros<mat>(K*T,K*T);
  for(int s=0; s<T; s++){
    ts_G_mat_prod_A_0_to_tau(span(K*s,K*(s+1)-1),span(K*s,K*T-1)) = GmatprodAstotau( K, vec_alpha_ij, tau, s, DeltaT, modA_mat);
  }
  return(ts_G_mat_prod_A_0_to_tau);
}

//===========================================================================================================
/* ***********************************
Function matHit
Description:
matHit a function that construct the observation matrix H_i(t) at time t.
Dependances: none
*/

arma::mat matHit(arma::vec X_i_t){
  // X_i_t : vector of observation
  int K = (int)X_i_t.size();
  int l=0; //loop variable
  int p=0; //loop variable

  for( int i=0; i< K; i++){
    if(!isnan(X_i_t[i])){
      p++;
    }
  }
  mat H_i_t = zeros(p,K); //initialisation of matrix H_i(t)
  // filling ofa matrix H_i(t)
  for( int i=0; i< K; i++){
    if(!isnan(X_i_t(i))){
      H_i_t(l,i) =1;
      l++;
    }
  }
  if(p==0){
    mat H_i_t = zeros(K,K);
  }
  return(H_i_t);
}

//===========================================================================================
//' After vectorising the vzector Yi, this function returns
//' a vector indicating missing values : 1 = observed value, 0 = missing value
//' 
//' @param Yi a matrix with possibly NAs
//' 
//' @return a vector of elements 0,1
//' @export
//' 
// [[Rcpp::export]]
arma::vec compoYiNA(arma::mat& Yi){
  int K = Yi.n_cols;
  int ni = Yi.n_rows;
  vec compo_Yi_NA = ones(K*ni);
  int p=0;
  for( int i=0; i< ni; i++){
    for( int j=0; j< K; j++){
      if(isnan(Yi(i,j))){
        compo_Yi_NA(p) =0;
      }
      p++;
    }
  }
  return(compo_Yi_NA);
}

//===========================================================================================
//' Function that returns Yi (a vector) without NAs values
//' 
//' @param Yi a vector with possibly NAs
//' 
//' @return a vector 
//' @export
//' 
// [[Rcpp::export]]
arma::vec YiwoNA(arma::vec Yi){
  int Ti = Yi.size();
  mat NAs = Yi.elem(find_nonfinite(Yi)); // repere NAs values in the vector Yi
  vec Yi_wo_NA = zeros(Ti-NAs.n_rows*NAs.n_cols);
  int p=0;
  for( int i=0; i< Ti; i++){
    if(!isnan(Yi(i))){
      Yi_wo_NA(p) =Yi(i);
      p++;
    }

  }
  return(Yi_wo_NA);
}
//===========================================================================================================

/* ***********************************
Function matNui
Description:
matNui a function that construct the matrix nu_t_j, the expectation of processes at time t_j
Dependance: none
*/

//===========================================================================================
//'  Function that constructs the matrix nu_t_j, the expectation of the processes at time t_j
//'  
//' @param nD an integer indicating the number of processes
//' @param tau_i vector of integers indicating times 
//' @param DeltaT double that indicates the discretization step
//' @param x0i model.matrix for baseline's submodel
//' @param xi model.matrix for change's submodel
//' @param alpha_mu0 a vector of parameters associated to the model.matrix for the baseline's submodel
//' @param alpha_mu a vector of parameters associated to the model.matrix for the change's submodel
//' @param G_mat_A_0_to_tau_i matrix containing  Prod(A_t)t=0,tau_i where A_t is the transition
//' matric containing at time t
//' 
//' @return a matrix
//' @export
//' 
// [[Rcpp::export]]
arma::mat matNui(int nD, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0,
                 arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i){
  // matNu_i : matrix of size xi.n_rows*nD containing expectation of processes
  int n_cols_xi = xi.n_cols;
  int mi=tau_i.size(); // number of observation
  int T = max(tau_i)+1;
  mat matNu_i = zeros(mi,nD);
  mat Mu_t = zeros(nD,1);
  int i = 0;
  for(int t=0; t< T; t++){
    if(t==0){
      Mu_t = x0i*alpha_mu0;
    }
    else{
      Mu_t = DeltaT*xi(span(t*nD,(t+1)*nD-1), span(0,n_cols_xi-1))*alpha_mu
      + G_mat_A_0_to_tau_i(span(0,nD-1),span(nD*(t-1),nD*(t-1)+nD-1))*Mu_t;
    }
    if(t ==tau_i(i)){
      matNu_i(span(i,i), span(0,nD-1)) = Mu_t.t();
      i++;
    }
  }
  return (matNu_i);
}
//===========================================================================================================
/* ***********************************
Function f_Yi_r_NA_by0
Description:
f_Yi_r_NA_by0 a function that replace NA by 0.0. Just for computatinal need
Dependences : none
*/

//===========================================================================================
//' Function that replaces NAs by 0.0 just for computatinal need
//' 
//' @param Yi a matrix with possibly NAs
//' 
//' @return a matrix 
//' @export
//' 
// [[Rcpp::export]]
arma::mat f_Yi_r_NA_by0(arma::mat& Yi){
  int K = Yi.n_cols;
  int ni = Yi.n_rows;
  vec compo_Yi_NA = ones(K*ni);
  int p=0;
  for( int i=0; i< ni; i++){
    for( int j=0; j< K; j++){
      if(isnan(Yi(i,j))){
        Yi(i,j) = 0.0;
      }
      p++;
    }
  }
  return(Yi);
}

//===========================================================================================================
/* ***********************************
Function YiNui
Description:
YiNui a function that compute the difference est une fonction (mat_Yi - mat_Nu_i)
delate missing values (NA) and return a vector
Dependances: matNui
*/

//===========================================================================================
//' Function that computes the difference (mat_Yi - mat_Nu_i), delates missing values (NAs) and 
//' returns a vector. mat_Yi is the outcomes and mat_Nu_i is the expectation
//'  
//' @param nD an integer indicating the number of processes
//' @param matrixP a matrix that matches markers to latent processes
//' @param tau a vector of integers indicating times 
//' @param tau_i a vector of integers indicating times for individual i
//' @param DeltaT double that indicates the discretization step
//' @param Yi a matrix of the outcomes
//' @param x0i model.matrix for baseline's submodel
//' @param xi model.matrix for change's submodel
//' @param alpha_mu0 a vector of parameters associated to the model.matrix for the baseline's submodel
//' @param alpha_mu a vector of parameters associated to the model.matrix for the change's submodel
//' @param G_mat_A_0_to_tau_i matrix containing  Prod(A_t)t=0,tau_i where A_t is the transition
//' matric containing at time t
//' 
//' @return a vector
//' @export
//' 
// [[Rcpp::export]]
arma::vec YiNui(int nD, arma::mat matrixP, arma::vec& tau, arma::vec& tau_i, double DeltaT, arma::mat& Yi, arma::mat& x0i, arma::colvec& alpha_mu0,
                arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i){
  // Yi : matrice of observation of subject i
  mat Nu_cp = matNui(nD, tau, DeltaT, x0i, alpha_mu0, xi, alpha_mu, G_mat_A_0_to_tau_i);
  mat Nu_cp_i = zeros(Yi.n_rows,nD);
  for(int i=0; i<(int)tau_i.size(); i++){
    Nu_cp_i.row(i) = Nu_cp.row(tau_i(i));
  }
  // Nu_cp : Expectation of overall times
  mat M =  Yi - Nu_cp_i*matrixP.t();
  return(YiwoNA(vectorise(M)));
}
//===========================================================================================================

// /*============================================================
// generate one random gaussian variable with mean and stddev
// ==============================================================*/
// double rand_normal(double mean, double stddev)
// {//Box muller method
//   static double n2 = 0.0;
//   static int n2_cached = 0;
//   if (!n2_cached)
//   {
//     double x, y, r;
//     do
//     {
//       x = 2.0*rand()/RAND_MAX - 1;
//       y = 2.0*rand()/RAND_MAX - 1;
//
//       r = x*x + y*y;
//     }
//     while (r == 0.0 || r > 1.0);
//     {
//       double d = sqrt(-2.0*log(r)/r);
//       double n1 = x*d;
//       n2 = y*d;
//       double result = n1*stddev + mean;
//       n2_cached = 1;
//       return result;
//     }
//   }
//   else
//   {
//     n2_cached = 0;
//     return n2*stddev + mean;
//   }
// }
//
// /*============================================================
// Generate multivariate gaussian vector with mean m and SD = LLt
// ==============================================================*/
// // [[Rcpp::export]]
// arma::mat mvnorm(int seed, arma::vec& m, arma::mat& SD){
//   mat L = chol(SD).t();
//   // mat L = SD.t();
//   int d = m.size();
//   vec xx = zeros(d); // for indpt standard gaussian vector
//   for( int i=0; i<d; i++){
//     srand (seed*33*i);
//     xx[i] = rand_normal(0.0,1.0);
//   }
//   return(m + L*xx);
// }

/*============================================================
 Generate multivariate gaussian vector with mean m and SD = LLt
==============================================================*/
arma::mat mvnorm(int seed, arma::vec m, arma::mat SD){
  Rcpp::Environment base("package:CInLPN");
  Rcpp::Function g = base["f_mvrnorm"];
  vec x = as<arma::vec>(wrap(g(Rcpp::_["seed"] = seed, Rcpp::_["m"] = m, Rcpp::_["sd"] = SD)));
  return(x);
}

// /*================================================================
// Function MC : Compute the prediction in real scale using Monte Carlo approch
//  Dependences : none
// ==================================================================*/
// // [[Rcpp::export]]
// arma::vec MC(int K, int nr, arma::vec& mu, arma:: mat& SD, List& knots, arma::vec& ParaTransformY, int degree){
//   Rcpp::Environment base("package:CInLPN");
//   Rcpp::Function f = base["R_MC"];
//   arma::vec yi = as<arma::vec>(wrap(f(K, nr, mu, SD, knots, ParaTransformY, degree)));
//   return(yi);
// }


//===========================================================================================
//' Function that transforms a vector to a matrix
//' 
//' @param y a vector 
//' @param K an integer indicating the number of columns of the returned matrix
//' @param m_i an integer indicating the number of rows of the returned matrix
//' 
//' @return a matrix 
//' @export
//' 
// [[Rcpp::export]]
arma::mat VecToMat(arma::vec& y, int K, int m_i){

  mat mat_y = zeros(m_i,K);
  for(int p = 0; p < m_i; p++){
    mat_y.row(p) = y(span(p*K, ((p+1)*K-1))).t();
  }

  return(mat_y);
}