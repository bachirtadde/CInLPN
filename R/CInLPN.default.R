# default function of CInLPN
CInLPN.default <- function(fixed_X0.models, fixed_DeltaX.models, randoms_X0.models, randoms_DeltaX.models, mod_trans.model, 
                           DeltaT, outcomes, nD, mapping.to.LP, link, knots=NULL, subject, data, Time,
                           makepred, MCnr,
                           paras.ini= NULL, indexparaFixeUser, paraFixeUser, maxiter, univarmaxiter, nproc = 1, 
                           epsa =0.00001, epsb = 0.00001, epsd= 0.001, print.info = FALSE,...)
{
  cl <- match.call()
  
  
  ################### created formated data ##########################
  data_F <- DataFormat(data=data, subject = subject, fixed_X0.models = fixed_X0.models,
                       randoms_X0.models = randoms_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                       randoms_DeltaX.models = randoms_DeltaX.models, mod_trans.model = mod_trans.model, 
                       outcomes = outcomes, nD = nD, link=link, knots = knots, 
                       Time = Time, DeltaT=DeltaT)
  
  
  K <- data_F$K #  number of markers
  vec_ncol_x0n <- data_F$vec_ncol_x0n # number of parameters on initial level of processes
  n_col_x <- ncol(data_F$x) # number of parameters on processes slope
  nb_RE <- data_F$nb_RE # number of random effects on the processes slope
  L <- ncol(data_F$modA_mat)
  ncolMod.MatrixY <- ncol(data_F$Mod.MatrixY)
  
  # ### creation of arguments:  Initialising parameters
  if(K>1 & is.null(paras.ini)){
    paras.ini <- f_paras.ini(data = data, outcomes = outcomes, mapped.to.LP = mapping.to.LP, fixed_X0.models = fixed_X0.models, fixed_DeltaX.models = fixed_DeltaX.models,  
                             randoms_DeltaX.models = randoms_DeltaX.models, nb_RE = nb_RE, mod_trans.model = mod_trans.model, 
                             subject = subject, Time = Time, link = link, knots = knots,
                             DeltaT = DeltaT, maxiter = univarmaxiter, epsd = epsd, nproc = nproc, print.info = print.info)
  }
  paras <- Parametre(K=K, nD = nD, vec_ncol_x0n, n_col_x, nb_RE, indexparaFixeUser = indexparaFixeUser, 
                     paraFixeUser = paraFixeUser, L = L, ncolMod.MatrixY = ncolMod.MatrixY, paras.ini=paras.ini)
  
  if_link <- rep(0,K)
  for(k in 1:K){
    if(link[k] !="linear") if_link[k] <- 1
  }
  
  # estimation
  est <- CInLPN.estim(K= K, nD = nD, mapping.to.LP = mapping.to.LP, data = data_F, if_link = if_link, DeltaT = DeltaT, 
                      paras = paras, maxiter = maxiter, nproc = nproc, epsa = epsa, epsb = epsb,
                      epsd = epsd, print.info = print.info)
  
  res <- list(conv = est$istop, v = est$v, best = est$b, ca = est$ca, cb = est$cb, rdm = est$rdm, 
              niter = est$iter, coefficients = est$coefficients, posfix = est$posfix)
  
  #   # fitted value and correlation matrix
  #   m <- length(data$tau)
  res$fitted.values <- NULL
  if(makepred & res$conv==1){
    ptm2 <- proc.time() 
    cat("Be patient, the program is computed the predictions for the observed marker(s) \n")
    cat("based on Monte Carlo integration and Newton Raphson algorithm for search of inverse of a transformation function \n ")
    cat(paste("Number on replique for MC integration : N=",MCnr),"\n")
    cat("Convergence criteria for the seach of inverse : eps=1.e-9 \n")
    cat("Execution may take a little time because the prediction calculation program is not yet parallelized. Thank you \n")
    #==================== prédictions =============================
    res$Marginal_Predict <- data_F$id_and_Time
    res$SubjectSpecific_Predict <- data_F$id_and_Time
    col <- colnames(res$Marginal_Predict)
    # colSS <- colnames(res$SubjectSpecific_Predict)
    if(requireNamespace("splines2", quietly = TRUE)){
    Predict <- pred(K = K, nD = nD, mapping = mapping.to.LP, paras = res$coefficients,
                    m_is= data_F$m_i, Mod_MatrixY = data_F$Mod.MatrixY, df= data_F$df,
                    x = data_F$x, z = data_F$z, q = data_F$q, nb_paraD = data_F$nb_paraD, x0 = data_F$x0, z0 = data_F$z0,
                    q0 = data_F$q0, if_link = if_link, tau = data_F$tau,
                    tau_is=data_F$tau_is, modA_mat = data_F$modA_mat, DeltaT=DeltaT, 
                    MCnr = MCnr, data_F$minY, data_F$maxY, data_F$knots, data_F$degree, epsPred = 1.e-9)
    }else{
      stop("Need package MASS to work, Please install it.")
    }
    kk <- 1
    for(k in 1: K){
      res$Marginal_Predict <- cbind(res$Marginal_Predict,data_F$Y[,k],Predict[,kk], 
                                    (data_F$Y[,k]-Predict[,kk]),Predict[,(kk+1):(kk+3)])
      
      res$SubjectSpecific_Predict <- cbind(res$SubjectSpecific_Predict,data_F$Y[,k],Predict[,(kk+4)], 
                                           (data_F$Y[,k]-Predict[,(kk+4)]),Predict[,c((kk+1),(kk+5),(kk+6))])
      
      col <- c(col,outcomes[k],paste(outcomes[k], "Pred", sep="."), paste(outcomes[k], "Res", sep="."),
               paste(outcomes[k], "tr", sep="-"), paste(outcomes[k], "tr.Pred", sep="-"),
               paste(outcomes[k], "tr.Res", sep="-"))
      kk <- kk+7
    }
    colnames(res$Marginal_Predict) <- col
    colnames(res$SubjectSpecific_Predict) <- col
    
    p.time2 <- proc.time() - ptm2
    cat("Prediction computation take :", p.time2[1], "\n")
  }
  
  ## Compute number of parameters  per componante of the processes network
  i1 <- 0
  res$length_para_mu0 <- length(which(res$posfix[(i1+1):(i1+ncol(data_F$x0))]==0))
  i1 <- i1 + ncol(data_F$x0)
  res$length_para_mu <- length(which(res$posfix[(i1+1):(i1+ncol(data_F$x))]==0))
  i1 <- i1 + ncol(data_F$x)
  res$length_para_RE <- length(which(res$posfix[(i1+1):(i1+data_F$nb_paraD)]==0))
  res$length_para_trY = ncol(data_F$Mod.MatrixY)
  
  ###names of estimates parameters 
  #colname transition matrix
  a <- NULL
  
  for(i in 1:nD){
    a <- c(a, paste("a",paste(i,1:nD,sep=""), sep="_"))
  }
  Col.matA <- NULL
  colmodA_mat <- colnames(data_F$modA_mat)
  for(i in 1: length(a)){
    Col.matA <- c(Col.matA, paste(a[i],colmodA_mat, sep="."))
  }
  # colname RE
  L <- NULL
  #   for(i in 1:K){
  npRE <- data_F$nb_paraD
  L <- paste("Chol.", 1:npRE,sep="")
  
  # colname  measurement error
  sig <- paste("sigma",outcomes, sep=".")
  # Measurement model parameters
  nb_para_trY <- ncol(data_F$Mod.MatrixY)
  ParaTransformY <- colnames((data_F$Mod.MatrixY))
  
  ##
  ##output related to the statistical model
  res$ns <- data_F$nb_subject
  res$N <- data_F$nb_obs
  res$nD <- data_F$nD
  res$K <- data_F$K
  # est$best : only estimates of non constraint parameters
  res$linkstype <- link
  res$linknodes <- data_F$knots
  res$nb_paraD <- data_F$nb_paraD
  res$tau = data_F$tau
  res$outcomes <- outcomes
  res$covariates <- data_F$all.pred
  res$Time <- Time 
  res$DeltaT <- DeltaT
  res$colnames <- c(colnames(data_F$x0), colnames(data_F$x), L, Col.matA, sig,ParaTransformY)
  res$coefficients <- as.matrix(res$coefficients)
  rownames(res$coefficients) <- res$colnames
  colnames(res$coefficients) <- "Coef."
  
  
  res$loglik <- est$fn.value
  res$AIC <- -2*est$fn.value + 2*length(est$b)
  res$BIC <- -2*est$fn.value + log(res$N)*length(est$b)
  
  
  ##output related to the iteration process
  # res$istop: binary indicator of if convergence is reached
  # res$iter: indicates the number of iterations
  # res$ca : convergence criteria related to the stability of the parameters
  # res$cb : convergence criteria related to the stability of the likelihood
  # res$rdm : convergence criteria related to the inversibility of the Hessian matrix
  
  
  
  ##output related big data like predictions, stocked matrix
  # res$fitted.values
  res$modA_mat = data_F$modA_mat
  
  res$call <- match.call()
  class(res) <- 'CInLPN'
  res
}
