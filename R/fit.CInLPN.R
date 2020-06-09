#' Marginal and subject-specific predictions for CInLPN objects
#'
#' @param object CInLPN object
#' @param newdata dataset
#' @param MCnr an integer that gives the number of Monte Carlo iterations
#' @param \dots optional parameters
#'
#' @return list of marginal and subject-specific predictions
#' @export
fit.CInLPN <- function(object, newdata, MCnr = 10, ...){
  model <- object
  cl <- match.call()
  if(missing(model)) stop("The argument model should be specified")
  if(class(model)!="CInLPN") stop("argument model must be a CInLPN object")
  x <- model$call
  if(missing(newdata)) stop("The argument newdata should be specified")
  
  ### identification of all components and sub-models of the model
  
  ## components of structural model
  fixed_X0 <- as.formula(x$structural.model$fixed.LP0)
  fixed_DeltaX <-  as.formula(x$structural.model$fixed.DeltaLP)
  randoms_DeltaX <-  as.formula(x$structural.model$random.DeltaLP)
  mod_trans <-  as.formula(x$structural.model$trans.matrix)
  DeltaT <- model$DeltaT
  
  # components of measurement model
  link <- as.formula(x$measurement.model$link.functions$links)
  knots <- as.formula(x$measurement.model$link.functions$knots)
  
  ## subject, Time
  subject <- x$subject
  Time <- x$Time
  
  
  ### checking newdata format
  
  colnames<-colnames(newdata)
  # if(missing(DeltaT) || DeltaT < 0 ) stop("DeltaT of the model must not be  null or negative")
  if(!(subject%in%colnames))stop("Subject should be in the data")
  if(!(Time %in% colnames)) stop("time should be in the data")
  if(!all(round((newdata[,Time]/DeltaT)-round(newdata[,Time]/DeltaT),8)==0.0))stop(paste("time must be a multiple of", DeltaT, sep = " "))
  if(dim(unique(newdata))[1] != dim(newdata)[1]) stop("Some rows are the same in the dataset, perhaps because of a too large discretization step")
  
  ### pre-processing of data
  
  ### outcomes and latent processes ####
  outcome <- as.character(attr(terms(fixed_DeltaX),"variables"))[2]
  outcomes_by_LP<-strsplit(outcome,"[|]")[[1]]
  nD <- length(outcomes_by_LP) # nD: number of latent process
  
  outcomes <- NULL
  mapping.to.LP <- NULL
  for(n in 1:nD){
    outcomes_n <- strsplit(outcomes_by_LP[n],"[+]")[[1]]
    outcomes_n <-as.character(sapply(outcomes_n,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE))
    outcomes_n <- unique(outcomes_n)
    outcomes <- c(outcomes, outcomes_n)
    mapping.to.LP <- c(mapping.to.LP, rep(n,length(outcomes_n)))
  }
  K <- length(outcomes)
  all.Y<-seq(1,K)
  
  ### pre-processing of fixed effect 
  fixed_X0.models =strsplit(gsub("[[:space:]]","",as.character(fixed_X0)),"~")[[2]]
  fixed_X0.models<- as.vector(strsplit(fixed_X0.models,"[|]")[[1]]) 
  #
  fixed_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(fixed_DeltaX)),"~")[[3]]
  fixed_DeltaX.models<-strsplit(fixed_DeltaX.model,"[|]")[[1]]# chaque model d'effet fixe mais en vu de connaitre tous les pred.fixed du modele multi
  
  ### pre-processing of random effect
  randoms_X0.models <- rep("1",nD)
  randoms_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(randoms_DeltaX)),"~")[[2]]
  randoms_DeltaX.models<-strsplit(randoms_DeltaX.model,"[|]")[[1]]    
  
  #### pre-processing of  mod_trans transition matrix 
  mod_trans.model=strsplit(gsub("[[:space:]]","",as.character(mod_trans)),"~")[[2]]
  
  data_F <- DataFormat(data=newdata, subject = subject, fixed_X0.models = fixed_X0.models,
                       randoms_X0.models = randoms_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                       randoms_DeltaX.models = randoms_DeltaX.models, mod_trans.model = mod_trans.model, 
                       outcomes = outcomes, nD = nD, link=link, knots = knots, 
                       Time = Time, DeltaT=DeltaT)
  
  ### calling C++ function pred to compute fitted vallues of the outcomes from newdata
  if_link <- rep(0,K)
  for(k in 1:K){
    if(link[k] !="linear") if_link[k] <- 1
  }
  
  Marginal_Predict <- data_F$id_and_Time
  SubjectSpecific_Predict <- data_F$id_and_Time
  col <- colnames(Marginal_Predict)
  if(requireNamespace("splines2", quietly = TRUE)){
    Predict <- fit(K = K, nD = nD, mapping = mapping.to.LP, paras = model$coefficients,
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
    Marginal_Predict <- cbind(Marginal_Predict,data_F$Y[,k],Predict[,kk], 
                              (data_F$Y[,k]-Predict[,kk]),Predict[,(kk+1):(kk+3)])
    
    SubjectSpecific_Predict <- cbind(SubjectSpecific_Predict,data_F$Y[,k],Predict[,(kk+4)], 
                                     (data_F$Y[,k]-Predict[,(kk+4)]),Predict[,c((kk+1),(kk+5),(kk+6))])
    
    col <- c(col,outcomes[k],paste(outcomes[k], "Pred", sep="."), paste(outcomes[k], "Res", sep="."),
             paste(outcomes[k], "tr", sep="-"), paste(outcomes[k], "tr.Pred", sep="-"),
             paste(outcomes[k], "tr.Res", sep="-"))
    kk <- kk+7
  }
  colnames(Marginal_Predict) <- col
  colnames(SubjectSpecific_Predict) <- col
  
  ### returning fitted values
  return(list(Marginal_Predict = Marginal_Predict, SubjectSpecific_Predict = SubjectSpecific_Predict))
}