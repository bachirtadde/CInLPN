#step 1  root computation of estimation
#' Title
#'
#' @param K number of outcomes
#' @param nD number of latent processes
#' @param mapping.to.LP indicates which outcome measured which latent process, it is a mapping table between
#'  outcomes and latents processes
#' @param data indicates the data frame containing all the variables for estimating the model
#' @param if_link indicates if non linear link is used to transform an outcome
#' @param DeltaT indicates the discretization step
#' @param paras initial values for parameters
#' @param maxiter maximum iteration
#' @param nproc number of processor to be used for running this package, default value is 1
#' @param epsa threshold for the convergence criterion on the parameters, default value is 1.e-4
#' @param epsb threshold for the convergence criterion on the likelihood, default value is 1.e-4
#' @param epsd threshold for the convergence criterion on the derivatives, default value is 1.e-3
#' @param print.info to print information during the liklihood optimization, default value is FALSE
#'
#' @return CInLPN object

CInLPN.estim <- function(K, nD, mapping.to.LP, data, if_link = if_link, DeltaT=1.0, paras, 
                         maxiter = 500, nproc = 1, epsa =0.0001, epsb = 0.0001,epsd= 0.001, print.info = FALSE){
  cl <- match.call()
  #  non parall Optimisation 
  # package loading
  # if(requireNamespace("marqLevAlgParallel", quietly = TRUE)){
    temp <- try(marqLevAlg(b = paras$paraOpt, fn = Loglik, nproc = nproc, .packages = NULL, epsa=epsa, epsb=epsb, epsd=epsd,
                           maxiter=maxiter, print.info = print.info, minimize = FALSE,
                           DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
                           K = K, nD = nD, mapping =  mapping.to.LP, m_is = data$m_i, if_link = if_link,
                           Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
                           x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
                           x0 = data$x0, z0 = data$z0, q0 = data$q0,tau = data$tau, tau_is=data$tau_is,
                           modA_mat = data$modA_mat)
    ,silent = FALSE)
    if(inherits(temp,'try-error')){
      est <- list(istop=20, v=rep(0,length=((length(paras$paraOpt))*((length(paras$paraOpt)+1)/2))) ,
                  fn.value=100000000, b=paras$paraOpt, ca=1,cb=1,rdm=1,ier=-1)
    }else{
      est <- temp
    }
  # }else{
  #   stop("Package MarqLevAlgParallel required for the optimization process")
  #   
  # }
  
  # (res <- Loglik(paraOpt = paras$paraOpt, DeltaT=DeltaT, paraFixe = paras$paraFixe, posfix = paras$posfix,
  #                K = K, nD = nD, mapping = mapping.to.LP, m_is = data$m_i, if_link = if_link,
  #                Mod_MatrixY = data$Mod.MatrixY, Mod_MatrixYprim = data$Mod.MatrixYprim, df=data$df,
  #                x = data$x, z = data$z, q = data$q, nb_paraD = data$nb_paraD,
  #                x0 = data$x0, z0 = data$z0, q0 = data$q0,tau = data$tau, tau_is=data$tau_is,
  #                modA_mat = data$modA_mat))
  
  #   (res <- pred(paraOpt = paras$paraOpt,paraFixe = paras$paraFixe, 
  #                  posfix = paras$posfix, DeltaT=DeltaT, K = K, nD = nD, mapping = mapping.to.LP, 
  #                  modA_mat = data$modA_mat, m_is = data$m_i, Y = data$Y, 
  #                  x = data$x, z = data$z, q = data$q, x0 = data$x0, 
  #                  z0 = data$z0, q0 = data$q0, nb_paraDu= data$nb_paraDu, 
  #                  nb_paraDw= data$nb_paraDw, tau = data$tau, tau_is=data$tau_is))    
  
  #estimating para + fixed para
  para <- paras$para
  para[which(paras$posfix==0)] <- est$b
  est$coefficients <- para
  est$posfix <- paras$posfix
  est
}
