#' Function to initialize parameters kin multivariate CInLPN model
#'
#' @param K number of the markers
#' @param nD number of the latent processes
#' @param vec_ncol_x0n vector of number of columns of model.matrix for baseline's submodel
#' @param n_col_x number of overall columns of model.matrix for change's submodel
#' @param nb_RE number of random effects
#' @param stochErr indicates if the structural model contain stochastique components
#' @param indexparaFixeUser position of parameters to be constrained
#' @param paraFixeUser values associated to the index of parameters to be constrained
#' @param L number of columns of model.matrix for temporal infuences model
#' @param paras.ini initial values for parameters, default values is NULL
#' @param ncolMod.MatrixY vector of number of columns of model.matrix for transformation submodel
#'
#' @return a list
#' 
#' 
Parametre <- function(K, nD, vec_ncol_x0n, n_col_x, nb_RE, stochErr=FALSE, indexparaFixeUser =NULL,
                      paraFixeUser=NULL, L = 1, paras.ini, ncolMod.MatrixY){
  cl <- match.call()
  #   require(MASS)
  #initialisation des parametres
  # L = number of parameters for each coefficient a of matrix A
  # K = number of outcomes
  #======================================================================================
  nb_paraD = nb_RE*(nb_RE+1)/2
  indexparaFixeForIden <- NULL
  # if user not specified initial parameters
  
  if(is.null(paras.ini)){
    p <- 0 # position in the initialize parameters
    cpt1 <- 0 # compteur pour tous les paras
    cpt2<-0 # compteur de boucle
    #alpha_mu0
    alpha_mu0 <-rep(.1,sum(vec_ncol_x0n))
    p <- p+ sum(vec_ncol_x0n)
    index_paraFixe_mu0_constraint <- NULL
    for(n in 1:nD){
      alpha_mu0[(cpt2+1)] <- 0
      cpt2 <- cpt2 + vec_ncol_x0n[n]
      cpt1 <- cpt1 + vec_ncol_x0n[n]
    }
    
    #alpha_mu
    alpha_mu <-rep(.3,n_col_x)
    p <- p+n_col_x
    cpt1 <- cpt1 + n_col_x
    
    
    alpha_D <-rep(.1,nb_paraD)
    to_nrow <- nb_RE
    i_alpha_D <- 0
    index_paraFixeDconstraint <- NULL
    for(n in 1:nD){
      alpha_D[i_alpha_D+1] <- 1
      i_alpha_D <- i_alpha_D + to_nrow
      cpt1 <- cpt1 + to_nrow
      to_nrow <- to_nrow -1
    }
    p <- p+nb_paraD
    
    # para of transition matrix vec_alpha_ij
    vec_alpha_ij <- rep(0.4, L*nD*nD)
    cpt1 <- cpt1 + L*nD*nD
    p <- p + L*nD*nD
    #paraB
    paraB <- NULL
    if(stochErr==TRUE){
      paraB <- rep(.15,nD)
      cpt1 <- cpt1 + nD
      p <- p + nD
    }
    #paraSig
    paraSig <- rep(.5, K)
    cpt1 <- cpt1 + K
    p <- p + K
    ### parameters of the link function
    ParaTransformY <- rep(1, ncolMod.MatrixY)
    
    cpt1 <- cpt1 + ncolMod.MatrixY
    p <- p + ncolMod.MatrixY
  }
  
  
  
  # if user specified initial parameters
  if(!is.null(paras.ini)){
    p <- 0 # position in the initialize parameters
    cpt1 <-0 # counter for parameterd
    cpt2<-0 # loop counter
    #alpha_mu0
    alpha_mu0 <- paras.ini[(p+1):(p+sum(vec_ncol_x0n))]
    p <- p+ sum(vec_ncol_x0n)
    index_paraFixe_mu0_constraint <-NULL
    for(n in 1:nD){
      alpha_mu0[(cpt2+1)] <- 0
      cpt2 <- cpt2 + vec_ncol_x0n[n]
      cpt1 <- cpt1 + vec_ncol_x0n[n]
    }
    paraFixe_mu0_constraint <- rep(1,nD)
    #alpha_mu
    alpha_mu <- paras.ini[(p+1):(p+n_col_x)]
    p <- p+n_col_x
    cpt1 <- cpt1 + n_col_x
    #alpha_D parameters for cholesky of all random effects
    alpha_D <- paras.ini[(p+1):(p+nb_paraD)]
    to_nrow <- nb_RE
    i_alpha_D <- 0
    index_paraFixeDconstraint <- NULL
    for(n in 1:nD){
      alpha_D[i_alpha_D+1] <- 1
      i_alpha_D <- i_alpha_D + to_nrow
      cpt1 <- cpt1 + to_nrow
      to_nrow <- to_nrow -1
    }
    p <- p+nb_paraD
    paraFixeDconstraint <- rep(1,nD)
    # para of transition matrix vec_alpha_ij
    vec_alpha_ij <- paras.ini[(p+1):(p + L*nD*nD)]
    p <- p + L*nD*nD
    cpt1 <- cpt1 + L*nD*nD
    # paraB
    paraB <- NULL
    if(stochErr==TRUE){
      
      paraB <- paras.ini[(p+1):(p + nD)]
      p <- p + nD
      cpt1 <- cpt1 + nD
    }
    #paraSig
    paraSig <- paras.ini[(p+1):(p + K)]
    p <- p + K
    cpt1 <- cpt1 + K
    ### para of link function
    ParaTransformY <- paras.ini[(p+1):(p + ncolMod.MatrixY)]
    cpt1 <- cpt1 + ncolMod.MatrixY
    p <- p + ncolMod.MatrixY
  }
  
  #final vector of initial parameters
  paras <- c(alpha_mu0, alpha_mu, alpha_D, vec_alpha_ij,  paraB, paraSig, ParaTransformY )
  
  #initialisation
  #   paraOpt <- paras
  posfix <- rep(0,length(paras)) # 0 = non fixe 1 = fixe # initialisation
  # constraining of parameters==============
  indexFixe <- indexparaFixeForIden
  if(!is.null(indexparaFixeUser)){
    if(length(indexparaFixeUser) != length(paraFixeUser)){
      stop("The length of paraFixe does not correspond with the length of indexparaFixe")
    }
    indexFixe <- sort(unique(c(indexFixe,indexparaFixeUser)))
  }
  paraFixe <- rep(NA, length(posfix))
  if(!is.null(paraFixeUser)){
    paraFixe[c(indexparaFixeUser)]<- paraFixeUser
  }
  paraFixe[index_paraFixe_mu0_constraint]<- rep(0,K)
  paraFixe[index_paraFixeDconstraint]<- rep(1,K)
  if(sum(!is.na(paraFixe))==0){
    paraFixe = -1
    paraOpt <- paras
  }else{
    paraFixe <- paraFixe[!is.na(paraFixe)]
    posfix[indexFixe] <- 1 # fixiation des paras d'indexes dans indexparaFixe
    paras[indexFixe] <- paraFixe
    paraOpt <- paras[-indexFixe]
  }
  
  return(list(para = paras, paraOpt = paraOpt, paraFixe = paraFixe, posfix = posfix, L = L))
}


#' Initialisation of parameters
#'
#' @param data indicates the data frame containing all the variables for estimating the model
#' @param outcomes names of the outcomes
#' @param mapped.to.LP indicates which outcome measured which latent process, it is a mapping table between
#'  outcomes and latents processes
#' @param fixed_X0.models fixed effects in the submodel for the baseline level of processes
#' @param fixed_DeltaX.models a two-sided linear formula object for specifying the response outcomes (one the left part of ~ symbol) 
#' and the covariates with fixed-effects (on the right part of ~ symbol) 
#' @param randoms_DeltaX.models random effects in the submodel for change over time of latent processes
#' @param randoms_X0.models random effects in the submodel for the baseline level of processes
#' @param nb_RE number of random effects
#' @param mod_trans.model model for elements of the temporal transition matrix, which captures 
#' the temporal influences between latent processes
#' @param subject indicates the name of the covariate representing the grouping structure
#' @param Time indicates the name of the covariate representing the time
#' @param link indicates link used to transform outcome
#' @param knots indicates position of knots used to transform outcomes 
#' @param DeltaT indicates the discretization step
#' @param maxiter maximum iteration
#' @param epsa threshold for the convergence criterion on the parameters, default value is 1.e-4
#' @param epsb threshold for the convergence criterion on the likelihood, default value is 1.e-4
#' @param epsd threshold for the convergence criterion on the derivatives, default value is 1.e-3
#' @param nproc number of processor to be used for running this package
#' @param print.info  to print information during the liklihood optimization, default value is FALSE
#'
#' @return a list

f_paras.ini <- function(data, outcomes, mapped.to.LP, fixed_X0.models, fixed_DeltaX.models, randoms_DeltaX.models, 
                        randoms_X0.models, nb_RE, mod_trans.model, subject, 
                        Time, link, knots, DeltaT, maxiter = 25, epsa = .0001, epsb = .0001,
                        epsd = .0001, nproc = 1, print.info = TRUE)
{
  cl <- match.call()
  
  K <- length(unlist(outcomes))
  nD <- length(fixed_X0.models)
  paras.ini <-  NULL
  para.fixed_X0 <- NULL
  para.fixed_DeltaX <- NULL
  para.RE <- NULL
  para.trans <- NULL
  para.Sig <- NULL
  ParaTransformY <- NULL
  paras.ini <- list()
  for(k in 1 : K){
    data <- data[!is.na(data[,outcomes[k]]),]
    fixed_DeltaX <- as.formula(paste(outcomes[k],"~",fixed_DeltaX.models[mapped.to.LP[k]], sep=""))
    fixed_X0 <- as.formula(paste("~", fixed_X0.models[mapped.to.LP[k]], sep=""))
    n_col_x_k <- ncol(model.matrix(object = fixed_DeltaX,data = data ))
    n_col_x0_k <- ncol(model.matrix(object = fixed_X0,data = data ))
    mod_randoms_DeltaX <- as.formula(paste("~", randoms_DeltaX.models[mapped.to.LP[k]], sep=""))
    mod_trans <- as.formula(paste("~", mod_trans.model, sep=" "))
    indexparaFixeUser <- c(1,(n_col_x0_k+n_col_x_k+1))
    paraFixeUser <- c(0,1)
    
    structural.model <- list(fixed.LP0 = fixed_X0,
                             fixed.DeltaLP = fixed_DeltaX,
                             random.DeltaLP = mod_randoms_DeltaX, 
                             trans.matrix = mod_trans, 
                             delta.time = DeltaT)
    measurement.model <- list(link.functions = list(links = link[k], knots = knots[k]))
    
    parameters = list(paras.ini = NULL, Fixed.para.index = indexparaFixeUser, Fixed.para.values = paraFixeUser)
    
    option = list(nproc = nproc, print.info = print.info, maxiter = maxiter)
    
    mod <- CInLPN(structural.model = structural.model, measurement.model = measurement.model, parameters = parameters,
                  option = option, Time = Time, subject = subject, data = data)
    L <- ncol(mod$modA_mat)
    
    ## compute number of paramter per component
    coefficients <- as.vector(mod$coefficients)
    i1 <- 0
    if(k==1 | (k>1 && (mapped.to.LP[k-1]!= mapped.to.LP[k]))){
      para.fixed_X0 <- c(para.fixed_X0, coefficients[(i1+1):(i1+n_col_x0_k)])
      i1 <- i1+n_col_x0_k
      para.fixed_DeltaX <- c(para.fixed_DeltaX, coefficients[(i1+1):(i1+n_col_x_k)])
      i1 <- i1 + n_col_x_k + mod$nb_paraD
    }
    else{
      i1 <- i1 + n_col_x0_k + mod$nb_paraD
    }
    
    if(k==1){
      para.trans <- c(para.trans, coefficients[(i1+1):(i1+L)])
    }
    if(k!=1 && (mapped.to.LP[k-1]!= mapped.to.LP[k]) ){ # 
      para.trans <- c(para.trans, rep(0,nD*L), coefficients[(i1+1):(i1+L)])
    }
    i1 <- i1+L
    
    para.Sig <- c(para.Sig,  mod$b[(i1+1):(i1+1)])
    i1 <- i1+1
    
    ParaTransformY <- c(ParaTransformY,  coefficients[(i1+1):(i1+mod$length_para_trY)])
    i1 <- i1+mod$length_para_trY
  }
  para.RE <- rep(0.1, (nb_RE*(nb_RE+1)/2))
  paras.ini <- c(para.fixed_X0, para.fixed_DeltaX, para.RE, para.trans, para.Sig, ParaTransformY) 
  return(paras.ini)
}

