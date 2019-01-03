arma::vec InnerProd(arma::vec v1, arma::vec v2);
arma::mat KmatDiag(arma::vec& Kvector);
arma::mat DparChol(int q, arma::vec& qvector);
double aijt(int t, arma::vec alpha_ijl,  arma::mat modA_mat);
arma::vec vecaijt( int K, int t, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat);
arma::mat ConstrA(int K, int t, double DeltaT, arma::vec vec_alpha_ij, arma::mat modA_mat);
arma::mat ProdA(int K, int t2, int t1, double DeltaT, arma::vec& vec_alpha_ij, arma::mat& modA_mat);
arma::mat GmatA0totaui( int K, arma::vec& vec_alpha_ij, arma::vec& tau_i, double DeltaT, arma::mat modA_mat);
arma::mat GmatprodAstotau( int K, arma::vec& vec_alpha_ij, arma::vec& tau, int t_ini, double DeltaT, arma::mat modA_mat);
arma::mat tsGmatprodA0totau(int K, arma::vec& vec_alpha_ij, arma::vec& tau, double DeltaT, arma::mat modA_mat);
arma::mat matHit(arma::vec X_i_t);
arma::vec compoYiNA(arma::mat& Yi);
arma::vec YiwoNA(arma::vec Yi);
arma::mat matNui(int nD, arma::vec& tau_i, double DeltaT, arma::mat& x0i, arma::colvec& alpha_mu0,
                 arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i);

arma::vec YiNui(int nD, arma::mat matrixP, arma::vec& tau, arma::vec& tau_i, double DeltaT, arma::mat& Yi, arma::mat& x0i, arma::colvec& alpha_mu0,
                arma::mat& xi, arma::colvec& alpha_mu, arma::mat& G_mat_A_0_to_tau_i);

arma::mat transformY(arma::mat& Y, arma::colvec& paraEtha0, arma::colvec& paraEtha1);
arma::vec vectorise(arma::mat& M);
arma::mat mvnorm(int seed, arma::vec m, arma::mat SD);
arma::mat VecToMat(arma::vec& y, int K, int m_i);
arma::mat f_inv_mat(arma::mat& B);
int f_mat_print( arma::mat& B);

