#' Causal Inference in a Latent Processes Network
#' 
#' This function estimates a dynamic model based on D latent processes observed through 
#' K longitudinal markers (K >= D) which includes temporal relationships between latent processes. 
#' The structural model for the D latent processes combines a
#' multivariate linear mixed model and system of difference equations to model trajectories of the processes and 
#' their temporal influences other time. The temporal influences between processes are captured through a time-dependent transition matrix that can be 
#' adjusted for covariates. The model includes as a special case, a standard multivariate mixed model in which associations between processes are not 
#' captured by temporal influences but by correlated random effects.
#' Longitudinal markers are related to their latent processes through equations of observation that involve parameterized 
#' link functions (parameters are estimated simultaneously with  those of the structural model). 
#' The link function can be linear for Gaussian marker or non linear (using a basis of I-splines) for
#' non Gaussian markers.
#' All parameters are estimated simultaneously in a maximum likelihood framework 
#' See the methodological paper available at :http://arxiv.org/abs/1806.03659
#' 
#' @param structural.model a list of 5 arguments used to specify the structural model: 
#' 
#' \code{structural.model$fixed.LP0}{ allows to specify the fixed effects in the submodel for the baseline level of processes.
#' Note that there is no need to specify a random effect model for the baseline level of processes as we
#' systematically set a process-specific random intercept (with variance fixed to 1 for identifiability purpose). 
#' For identifiability purposes, the mean intercepts are fixed to 0 (not estimated).}
#' 
#' \code{structural.model$fixed.DeltaLP}{ a two-sided linear formula object for specifying the response outcomes (one the left part of ~ symbol) 
#' and the covariates with fixed-effects (on the right part of ~ symbol) 
#' in the submodel for change over time of latent processes.}
#' 
#' \code{structural.model$random.DeltaLP}{ allow to specify the random effects in the submodel for change over time of latent processes.}
#' 
#' \code{structural.model$trans.matrix}{ allow to specify a model for elements of the transition matrix, which captures 
#' the temporal influences between latent processes.}
#' 
#' \code{structural.model$delta.time}{ indicates the discretisation step to be used for latent processes}
#' 
#' @param measurement.model is a list of arguments detailed below used to specify the measurement model: 
#' 
#' \code{measurement.model$link}{ indicates the link functions to be used to transform the outcomes. 
#' It takes values in "linear" for a linear transformation and "n-type-d" for a I-splines transformation 
#' where "n" indicates the number of nodes, "type" (which takes values in "quant", "manual", "equi") indicates 
#' where the nodes are placed, and "d" indicates the degree of the I-splines.}
#' 
#' \code{measurement.model$knots}{ argument indicates if necessary the place of knots (when placed manually with "manual"), 
#' default value is NULL}
#' 
#' @param parameters a list of 3 arguments about parameters of the models (e.g., initial paremeters, parameters one would like to fix, etc.)
#' 
#' \code{parameters$paras.ini}{ indicates initial values for parameters, default values is NULL.}
#' 
#' \code{parameters$Fixed.para.indix}{ indicates the positions of parameters to be constrained.}
#' 
#' \code{parameters$Fixed.para.values}{ indicates the values associated to the index of parameters to be constrained. }
#' 
#' @param option a list of 4 arguments for the optimization procedure:
#' 
#' \code{option$epsa}{ threshold for the convergence criterion on the parameters}
#' 
#' \code{option$epsb}{ threshold for the convergence criterion on the likelihood}
#' 
#' \code{option$epsc}{ threshold for the convergence criterion on the derivatives}
#' 
#' \code{option$MCnr}{ number of replicates  to compute the predictions in the real scales of the outcomes (backward transformation because of the link functions)}
#' @param Time indicates the name of the covariate representing the time 
#' @param subject indicates the name of the covariate representing the grouping structure
#' @param data indicates the data frame containing all the variables for estimating the model.
#' @param  \dots other optional arguments
#' @return ---
#' @export
#'
#' @examples
#'         
#' ### example 1
#' Delta <- 1
#' paras.ini <- c(0.000, 0.059, 0.000, 0.163, -0.050, -0.153, 1.000, 0.000, 0.322, 0.000, 1.000, 
#'                0.000, 0.077, 0.139, 0.000, 0.177, -0.354, 0.114, 0.116, -0.090, 0.287, 0.554,
#'                1.107, 1.889, 0.881, 1.329)
#' indexparaFixeUser <- c(1,3, 6+c(1, 2, 4, 5, 6, 9))
#' paraFixeUser <- c(0, 0, 1, 0, 0, 1, 0, 0)
#' mod1 <- CInLPN(structural.model = list(fixed.LP0 = ~ 1 + C2 | 1 + C2,
#'                                       fixed.DeltaLP = L2 | L3  ~ 1 | 1 ,
#'                                       random.DeltaLP = ~ 1|1,
#'                                       trans.matrix = ~ 1,
#'                                       delta.time = Delta),
#'               measurement.model = list(link.functions = list(links = c(NULL,NULL),
#'                                                              knots = list(NULL, NULL))),
#'               
#'               parameters = list(paras.ini = paras.ini, Fixed.para.index = indexparaFixeUser, 
#'                                 Fixed.para.values = paraFixeUser),
#'               option = list(nproc = 1, print.info = TRUE, mekepred = TRUE, MCnr = 10, 
#'                             univarmaxiter = 7, epsa = 1e-5, epsb = 1e-4, epsd = 1e-2),
#'               Time = "time",
#'               subject = "id",
#'               data = data
#' )
#' 
#' summary(mod1)
#' 
#'  ### example 2
#'  library(splines)
#'  Delta <- 0.5
#'  paras.ini <- c(0.000,0.065, 0.000, 0.168, -0.054, 0.000, -0.119, -0.009, 1.000, 0.000,
#'                 0.473, 0.000, 1.000, 0.000, 0.057, -0.182, 0.000, 0.174, -0.523, 0.000,
#'                 0.000, 0.000, 0.073, 0.083, 0.119,  0.106, 0.015,  0.157,  0.054, 0.087,
#'                 -0.079, 0.000, 0.000, 0.000, 8.876, 0.287, 0.546, 0.531, 0.256, 1.100,
#'                 1.891,  0.846, 1.345)
#' indexparaFixeUser <- c(1,3, 8+c(1, 2, 4, 5, 6, 9,10+c(2:4,14:16)))
#' paraFixeUser <- c(0, 0, 1, 0, 0, 1, 0, 0, rep(0,6))
#' mod2 <- CInLPN(structural.model = list(fixed.LP0 = ~ 1 + C2 | 1 + C2,
#'                                      fixed.DeltaLP = L1 + L2| L3  ~ 1 + time| 1 + time,
#'                                      random.DeltaLP = ~ 1|1,
#'                                      trans.matrix = ~ 1 + bs(x = time, knots =c(2), 
#'                                                              intercept = F, degree = 2),
#'                                      delta.time = Delta),
#'              measurement.model = list(link.functions = list(links = c(NULL,NULL, NULL),
#'                                                             knots = list(NULL, NULL, NULL))),
#'              
#'              parameters = list(paras.ini = paras.ini, Fixed.para.index = indexparaFixeUser, 
#'                                Fixed.para.values = paraFixeUser),
#'              option = list(nproc = 2, print.info = TRUE, mekepred = TRUE, MCnr = 10, 
#'                            univarmaxiter = 7, epsa = 1e-5, epsb = 1e-4, epsd = 1e-2),
#'              Time = "time",
#'              subject = "id",
#'              data = data
#' )
#'
#' summary(mod2)
#'
#'#' 
#' \dontrun{
#' ### example 3
#' Delta <- 1
#' paras.ini <- NULL
#' indexparaFixeUser <- c(1,4,10+c(1,2,4,5,6,9, 10+c(1:4)))
#' paraFixeUser <- c(0,0,1,0,0,1,0,0, rep(0,4))

#' mod <-  CInLPN(structural.model = list(fixed.LP0 = ~ 1 + C1 + C2|1 + C1 + C2,
#'                                    fixed.DeltaLP = L1 + L2 | L3 ~1 + time|1 + time,
#'                                    random.DeltaLP = ~ 1|1,
#'                                    trans.matrix = ~ 1,
#'                                    delta.time = Delta),
#'            measurement.model = list(link.functions = list(links = c("4-equi-2", "linear",
#'                                                                     "4-equi-2"),
#'                                                           knots = list(NULL, NULL, NULL))),
#'            parameters = list(paras.ini = paras.ini, Fixed.para.index = indexparaFixeUser, 
#'                              Fixed.para.values = paraFixeUser),
#'            option = list(nproc = 1, print.info = FALSE, mekepred = TRUE, MCnr = 10, 
#'                          univarmaxiter = 7, epsa = 1e-5, epsb = 1e-5, epsd = 1e-5),
#'            Time = "time",
#'            subject = "id",
#'            data = data
#'          )
#'          
#' ### example 4
#' library(splines)
#' Delta <- 0.5
#' paras.ini <- NULL
#' indexparaFixeUser <- c(1,3, 8+c(1, 2, 4, 5, 6, 9,10+c(2:4,14:16)))
#' paraFixeUser <- c(0, 0, 1, 0, 0, 1, 0, 0, rep(0,6))
#' res <- CInLPN(structural.model = list(fixed.LP0 = ~ 1 + C2 | 1 + C2,
#'                                  fixed.DeltaLP = L1 | L2  ~ 1 + time| 1 + time,
#'                                  random.DeltaLP = ~ 1|1,
#'                                  trans.matrix = ~ 1 + bs(x = time, knots =c(2), 
#'                                                          intercept = F, degree = 2),
#'                                  delta.time = Delta),
#'          measurement.model = list(link.functions = list(links = c(NULL,NULL),
#'                                                         knots = list(NULL, NULL))),
#'          
#'          parameters = list(paras.ini = paras.ini, Fixed.para.index = indexparaFixeUser, 
#'                            Fixed.para.values = paraFixeUser),
#'          option = list(nproc = 2, print.info = TRUE, mekepred = TRUE, MCnr = 10, 
#'                        univarmaxiter = 7, epsa = 1e-5, epsb = 1e-4, epsd = 1e-2),
#'          Time = "time",
#'          subject = "id",
#'          data = data
#'   )
#'}
#'       


CInLPN <- function(structural.model, measurement.model, parameters, 
                   option, Time, subject, data,...){
  cl <- match.call()
  ptm <- proc.time()  
  cat("Be patient, CInLPN is running ... \n")
  
  ### check if all component of the model specification are well filled ####
  if(missing(structural.model))stop("The argument structural.model must be specified")
  if(missing(measurement.model))stop("The argument measurement.model must be specified")
  if(missing(parameters))stop("The argument parameters must be specified")
  if(missing(subject))stop("The argument subject must be specified")
  if(missing(Time))stop("The argument time must be specified")
  if(missing(data))stop("The argument data must be specified")
  
  if(is.null(structural.model$fixed.DeltaLP))stop("The argument structural.model$fixed.DeltaLP must be specified")
  if(is.null(structural.model$random.DeltaLP))stop("The argument structural.model$random.DeltaLP must be specified")
  if(is.null(structural.model$trans.matrix))stop("The argument structural.model$trans.matrix must be specified")
  if(is.null(structural.model$delta.time)){
    structural.model$delta.time <- 1
  }
  if(is.null(measurement.model$link.functions)){
    links <- NULL
    knots <- NULL
    measurement.model$link.functions =list(links = links, knots = knots)
  }
  if(is.null(parameters$Fixed.para.index))stop("The argument parameters$Fixed.para.index cannot be NULL")
  if(is.null(parameters$Fixed.para.values))stop("The argument parameters$Fixed.para.values cannot be NULL")
  
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
  # if(missing(DeltaT) || DeltaT < 0 ) stop("The discretization step DeltaT cannot be null or negative")
  if(!(subject%in%colnames))stop("Subject should be in the dataset")
  if(!(Time %in% colnames)) stop("Time variable should be indicated and should be in the dataset")
  if(!all(round((data[,Time]/DeltaT)-round(data[,Time]/DeltaT),8)==0.0))stop(paste("Time must be multiple of", DeltaT, sep = " "))
  if(dim(unique(data))[1] != dim(data)[1]) stop("Some rows are the same in the dataset, perhaps because of a too large discretisation step")
  
  
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
    if(is.null(outcomes_n)) stop("at least one marker must be specified for each latent process" )
    outcomes <- c(outcomes, outcomes_n)
    mapping.to.LP <- c(mapping.to.LP, rep(n,length(outcomes_n)))
  }
  if(!all(outcomes%in% colnames)) stop("outcomes must be in the data")
  K <- length(outcomes)
  all.Y<-seq(1,K)
  
  fixed_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(fixed_DeltaX)),"~")[[3]]
  fixed_DeltaX.models<-strsplit(fixed_DeltaX.model,"[|]")[[1]]# chaque model d'effet fixe mais en vu de connaitre tous les pred.fixed du modele multi
  
  if(nD !=length(fixed_DeltaX.models)) stop("The number of models does not correspond to the number of latent processes")
  
  if(nD > K){
    stop("There are too many latent processes compared to the indicated number of markers")
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
    stop("The number of models for initial latent processes does not correspond with the number of latent processes")
  }
  if(nD!=length(fixed_DeltaX.models)){
    stop("The number of models for the change over time of latent processes does not correspond with the number of latent processes")
  }
  
  ### traitement of transformation models ##
  if(is.null(link)){
    link <- rep("linear",K)
  }
  else if(length(link)!=K) stop("The number transformation links must be equal to the number of markers")
  
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
  cat("The program took:", p.time[1], "\n")
  est
}

