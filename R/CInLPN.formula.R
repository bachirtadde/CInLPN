#' Causal Inference in a Latent Processes Network
#' 
#' This package  estimates trajectories of multivariate latent processes observed through K longitudinal markers (K ≥ D) and
#' relationships  between latent processes. 
#' A structural model for processes network is defined for latent processes. This structural model combines 
#' multivariate linear mixed model and system of difference equations to model trajectories of the processes and 
#' their temporal influences other time. The temporal influences are captured through a time-depend transition matrix that can be 
#' adjusted for covariates. Instead of temporal influences, the variance-covariance matrix of random effects on processes can be suitably used to
#' modeling  (residual) correlations between latent processes.
#' Longitudinal markers are related to their latent processes through link function, which parameters are estimated simultaneously 
#' with  those of structural model. The link function can be linear for Gaussien marker or non linear (I-splines basis) for
#' non Gaussian markers.
#' All parameters are  estimated simultaneously in a maximum likelihood framework 
#' See the methodological paper available at :http://arxiv.org/abs/1806.03659
#' 
#' @param structural.model a list of arguments used to specify the structural model
#' The structural model contains  five component to specify it.
#' @param structural.model$fixed.LP0 allow to specify a fixed effects model on the baseline level of processes
#' Note that there is no need to specify random effect model for the baseline level of processes. As we
#' set systematicaly an random intercept with variance 1 for identifiability purpose
#' @param structural.model$fixed.DeltaLP a two-sided linear formula object for specifying the fixed-effects in the linear mixed 
#' model for initial level of latent processes. The response outcome is on the left of ~ symbol and the covariates are separated by + on the 
#' right of the ~ symbol. Fo identifiability purposes, the intercept of model are fixed to 0 (not estimated).
#' @param structural.model$random.DeltaLP allow to specify a random effects model on the processes slopes
#' @param structural.model$trans.matrix allow to specify a model for element of the transition matrix, which capture 
#' temporal influences between latent processes.
#' @param structural.model$delta.time to indicates the discretisation step to be used for latent processes
#' @param measurement.model a list of arguments used to specify the measurement model
#' @param measurement.model$link to indicate link to be used to transform markers. Takes values in ("linear", "n-type-d")
#' where "n" indicates the number of nodes, "type" : takes values in ("quant", "manual", "equi") to indicates how nodes are placed, 
#' and finaly argument "d" indicates the degre of I-splilnes to be used to construct the link function
#' @param measurement.model$knots argument to indicate if necessary the place of kodes, default value is NULL
#' @param parameters a list of arguments about parameters of the models (e.g., initial paremeters, parameters one would like to fix, etc.)
#' @param parameters$paras.ini to indicates initial values for parameters, default values is NULL
#' @param parameters$Fixed.para.indix to indicates the position of parameters to be constraint. 
#' @param parameters$Fixed.para.values to indicates the values associates to the index of parameters to be constraint. 
#' @param option a list of arguments for optimization purpose
#' @param option$epsa convergence criteria
#' @param option$epsb convergence criteria
#' @param option$epsc convergence criteria
#' @param option$MCnr number of replik  to compute the predictions in real scales
#' @param Time name of the covariate representing the time 
#' @param subject name of the covariate representing the grouping structure
#' @param data data frame containing the variables named in strutural model, time, subject.
#'
#' @return 0
#' @export
#'
#' @examples
#' ### example 1
#' Delta <- 1
#' paras.ini <- NULL
#' indexparaFixeUser <- c(1,4,10+c(1,2,4,5,6,9, 10+c(1:4)))
#' paraFixeUser <- c(0,0,1,0,0,1,0,0, rep(0,4))

#' mod <-  CInLPN(structural.model = list(fixed.LP0 = ~ 1 + C1 + C2|1 + C1 + C2,
#'                                    fixed.DeltaLP = L1 + L2 | L3 ~1 + time|1 + time,
#'                                    random.DeltaLP = ~ 1|1,
#'                                    trans.matrix = ~ 1,
#'                                    delta.time = Delta),
#'            measurement.model = list(link.functions = list(links = c("4-equi-2", "linear", "4-equi-2"),
#'                                                           knots = list(NULL, NULL, NULL))),
#'            parameters = list(paras.ini = paras.ini, Fixed.para.index = indexparaFixeUser, Fixed.para.values = paraFixeUser),
#'            option = list(nproc = 1, print.info = F, mekepred = T, MCnr = 10, univarmaxiter = 7, epsa = 1e-5, epsb = 1e-5, epsd = 1e-5),
#'            Time = "time",
#'            subject = "id",
#'            data = data
#'          )
#'          
#' ### example 2
#' library(splines)
#' Delta <- 0.5
#' paras.ini <- NULL
#' indexparaFixeUser <- c(1,3, 8+c(1, 2, 4, 5, 6, 9,10+c(2:4,14:16)))
#' paraFixeUser <- c(0, 0, 1, 0, 0, 1, 0, 0, rep(0,6))
#' res <- CInLPN(structural.model = list(fixed.LP0 = ~ 1 + C2 | 1 + C2,
#'                                  fixed.DeltaLP = L1 | L2  ~ 1 + time| 1 + time ,
#'                                  random.DeltaLP = ~ 1|1,
#'                                  trans.matrix = ~ 1 + bs(x = time, knots =c(2), intercept = F, degree = 2),
#'                                  delta.time = Delta),
#'          measurement.model = list(link.functions = list(links = c(NULL,NULL),
#'                                                         knots = list(NULL, NULL))),
#'          
#'          parameters = list(paras.ini = paras.ini, Fixed.para.index = indexparaFixeUser, Fixed.para.values = paraFixeUser),
#'          option = list(nproc = 2, print.info = T, mekepred = T, MCnr = 10, univarmaxiter = 7, epsa = 1e-5, epsb = 1e-4, epsd = 1e-2),
#'          Time = "time",
#'          subject = "id",
#'          data = data
#'   )
#'     
CInLPN <- function(structural.model, measurement.model, parameters, 
                   option, Time, subject, data,...){
  cl <- match.call()
  ptm <- proc.time()  
  cat("Be patient, causal Model is running ... \n")
  
  ### check if all component of the model soecification are well filled ####
  if(missing(structural.model))stop("The argument structural.model must be specified in any model")
  if(missing(measurement.model))stop("The argument measurement.model must be specified in any model")
  if(missing(parameters))stop("The argument parameters must be specified in any model")
  if(missing(subject))stop("The argument subject must be specified in any model")
  if(missing(Time))stop("The argument time must be specified in any model")
  if(missing(data))stop("The argument data must be specified in any model")
  
  if(is.null(structural.model$fixed.DeltaLP))stop("The argument structural.model$fixed.DeltaLP must be specified in any model")
  if(is.null(structural.model$random.DeltaLP))stop("The argument structural.model$random.DeltaLP must be specified in any model")
  if(is.null(structural.model$trans.matrix))stop("The argument structural.model$trans.matrix must be specified in any model")
  if(is.null(structural.model$delta.time)){
    structural.model$delta.time <- 1
  }
  if(is.null(measurement.model$link.functions)){
    links <- NULL
    knots <- NULL
    measurement.model$link.functions =list(links = links, knots = knots)
  }
  if(is.null(parameters$Fixed.para.index))stop("The argument parameters$Fixed.para.index must not be NULL")
  if(is.null(parameters$Fixed.para.values))stop("The argument parameters$Fixed.para.values must not be NULL")
  
  if(is.null(option$makepred)){
    option$makepred <- TRUE
  }
  if(is.null(option$MCnr)){
    option$MCnr <- 30
  }
  
  # if(is.null(option$parallel)){
  #   option$parallel <- FALSE
  # }
  if(is.null(option$maxiter)){
    option$maxiter <- 500
  }
  if(is.null(option$univarmaxiter)){
    option$univarmaxiter <- 25
  }
  if(is.null(option$nproc)){
    option$nproc <- 1
  }
  if(is.null(option$print.info)){
    option$print.info <- FALSE
  }
  if(is.null(option$epsa)){
    option$epsa <- 0.005
  }
  if(is.null(option$epsb)){
    option$epsb <- 0.005
  }
  if(is.null(option$epsd)){
    option$epsd <- 0.05
  }
  
  ### identification of model components #####
  ## components of structural model
  fixed_X0 <- structural.model$fixed.LP0
  fixed_DeltaX <- structural.model$fixed.DeltaLP
  randoms_DeltaX <- structural.model$random.DeltaLP
  mod_trans <- structural.model$trans.matrix
  DeltaT <- structural.model$delta.time
  
  # components of measurement model
  link <- measurement.model$link.functions$links
  knots <- measurement.model$link.functions$knots
  ## components of parameters initialisation
  indexparaFixeUser <- parameters$Fixed.para.index
  paraFixeUser <- parameters$Fixed.para.values
  paras.ini <- parameters$paras.ini
  ## component of option
  makepred <- option$makepred
  MCnr <- option$MCnr
  #parallel <- option$parallel
  maxiter <- option$maxiter  
  univarmaxiter <- option$univarmaxiter
  nproc <- option$nproc
  epsa <- option$epsa
  epsb <- option$epsb
  epsd <- option$epsd
  print.info <- option$print.info
  
  colnames<-colnames(data)
  # if(missing(DeltaT) || DeltaT < 0 ) stop("DeltaT of the model must not be  null or negative")
  if(!(subject%in%colnames))stop("Subject should be in colnames")
  if(!(Time %in% colnames)) stop("time must be set and must be in colonum names of the dataset")
  if(!all(round((data[,Time]/DeltaT)-round(data[,Time]/DeltaT),8)==0.0))stop(paste("time must be multiple of", DeltaT, sep = " "))
  if(dim(unique(data))[1] != dim(data)[1]) stop("Some rows are the same in the dataset, Perhaps because of the discretisation step")
  
  
  #### fixed effects pre-traitement ####
  
  ### for DeltaLP
  # if(missing(fixed_DeltaX)) stop("The argument fixed_DeltaX must be specified in any model")
  if(class(fixed_DeltaX)!="formula") stop("The argument fixed_DeltaX must be a formula")
  
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
    if(is.null(outcomes_n)) stop("at least one marker must be specified for a latent process" )
    outcomes <- c(outcomes, outcomes_n)
    mapping.to.LP <- c(mapping.to.LP, rep(n,length(outcomes_n)))
  }
  if(!all(outcomes%in% colnames)) stop("outcomes of the model must be  colname of the data frame")
  K <- length(outcomes)
  all.Y<-seq(1,K)
  
  fixed_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(fixed_DeltaX)),"~")[[3]]
  fixed_DeltaX.models<-strsplit(fixed_DeltaX.model,"[|]")[[1]]# chaque model d'effet fixe mais en vu de connaitre tous les pred.fixed du modele multi
  
  if(nD !=length(fixed_DeltaX.models)) stop("There are less/more models than number of latent processes")
  
  if(nD > K){
    stop("There are more latent processes than markers")
  }
  
  ### pre-traitement of fixed effect on initial levels of processes
  if(is.null(fixed_X0)){
    fixed_X0<- ~1
    fixed_X0.models <- rep("1",nD)
  }
  if(class(fixed_X0)!="formula") stop("The argument fixed_X0 must be a formula")
  
  fixed_X0.models =strsplit(gsub("[[:space:]]","",as.character(fixed_X0)),"~")[[2]]
  fixed_X0.models<- as.vector(strsplit(fixed_X0.models,"[|]")[[1]]) 
  for(nd in 1:nD){
    if(fixed_X0.models[nd]=="~-1") fixed_X0.models[nd] <-"~1" # au moins l'intcpt
  }
  
  ### pre-traitement of random effect on processes  intercept and slope
  #### randoms effet on DeltaLP 
  randoms_X0.models <- rep("1",nD)
  #### randoms effet on DeltaX
  if(missing(randoms_DeltaX)){
    randoms_DeltaX<- ~1
    randoms_DeltaX.models <- rep("1",nD)
  }
  if(class(randoms_DeltaX)!="formula") stop("The argument random must be a formula")
  randoms_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(randoms_DeltaX)),"~")[[2]]
  randoms_DeltaX.models<-strsplit(randoms_DeltaX.model,"[|]")[[1]]    
  
  #### traitement of  mod_trans: transition matrix##### 
  if(missing(mod_trans)){
    mod_trans <- ~ 1 # constant transition matrix
  } 
  if(class(mod_trans)!="formula") stop("The argument mod_trans must be a formula")
  mod_trans.model=strsplit(gsub("[[:space:]]","",as.character(mod_trans)),"~")[[2]]
  
  if(nD!=length(fixed_X0.models)){
    stop("Diference between numbers of models for initial latent processes and number of latent processes")
  }
  if(nD!=length(fixed_DeltaX.models)){
    stop("Diference between numbers of models for the slope of latent processes and number of latent processes")
  }
  
  ### traitement of transformation models ##
  if(is.null(link)){
    link <- rep("linear",K)
  }
  else if(length(link)!=K) stop("number transformation links must be equal to the number of markers")
  
  #### call of CInLPN.default function to compute estimation and predictions
  est <- CInLPN.default(fixed_X0.models = fixed_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, randoms_X0.models = randoms_X0.models, 
                        randoms_DeltaX.models = randoms_DeltaX.models, mod_trans.model = mod_trans.model, DeltaT = DeltaT , outcomes = outcomes,
                        nD = nD, mapping.to.LP = mapping.to.LP, link = link, knots = knots, subject = subject, data = data, Time = Time, 
                        makepred = option$makepred, MCnr = option$MCnr, paras.ini= paras.ini, paraFixeUser = paraFixeUser, indexparaFixeUser = indexparaFixeUser,  
                        maxiter = maxiter, univarmaxiter = univarmaxiter, nproc = nproc, epsa = epsa, epsb = epsb, epsd = epsd, 
                        print.info = print.info)
  est$call <- match.call()
  est$formula <- list(fixed_X0.models=fixed_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                      randoms_X0.models=randoms_X0.models, randoms_DeltaX.models=randoms_DeltaX.models, 
                      mod_trans.model = mod_trans.model)
  est$mapping.to.LP <- mapping.to.LP
  
  p.time <- proc.time() - ptm
  cat("The program took :", p.time[1], "\n")
  est
}

