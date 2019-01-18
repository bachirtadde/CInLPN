# =======================================================
# function that return 1 if all element of a vector is NA
# =======================================================
is_na_vec <- function(vec){
  cl <- match.call()
  size <- length(t(vec))
  is_na_vec <- 0
  for(i in 1: size){
    is_na_vec[i] <- as.numeric(is.na(vec[i]))
  }
  return(prod(is_na_vec))
}

###### link function
f.link <- function(outcomes, Y,link=NULL, knots = NULL, na.action = 'na.pass'){
  cl <- match.call()
  current.na.action <- options('na.action')
  options(na.action = na.action)
  if(requireNamespace("splines2", quietly = TRUE)){
    col <- colnames(Y)
    K <- length(col)
    if(is.null(link)) link <- rep("linear",K)
    if(!is.null(link) && length(gsub("[[:space:]]","",link)) != K){
      stop("The number of transformation links does not correspond to the number of markers")
    }else{
      link <- gsub("[[:space:]]","",link)
    }
    degree <- rep(NA,K)
    if(is.null(knots)) knots <- vector("list", K)
    linkSpe <-list()
    Mod.MatrixY <- NULL
    Mod.MatrixYprim <- NULL
    minY <- rep(NA,K)
    maxY <- rep(NA,K)
    df <- NULL
    colnamesY <- NULL
    colnamesYPrim <- NULL
    for(k in 1:K){
      minY[k] <- min(Y[,col[k]], na.rm = T)
      maxY[k] <- max(Y[,col[k]], na.rm = T)
      if(link[k]=="linear"){
        linkSpe[[k]] <- "-"
        Imat <- model.matrix(as.formula(paste("~1+",col[k])), data = Y, na.action= na.action)
        Imat[which(is.na(Imat[,2])),1] <- NA
        colnamesY <- c(colnamesY, paste(outcomes[k],"linear", seq(1,ncol(Imat)), sep = "."))
        colnamesYPrim <- c(colnamesYPrim, paste(outcomes[k],"linear", seq(1,(ncol(Imat)-1)), sep = "."))
        
        Mod.MatrixY <- cbind(Mod.MatrixY, Imat)
        Mod.MatrixYprim <- cbind(Mod.MatrixYprim, Imat[,1])
        
        df <-c(df, ncol(Imat))
        degree[k] <- 0 # conventionnellement
      }
      else{
        linkSpe[[k]] <- strsplit(gsub("[[:space:]]","",link[k]),"[-]")[[1]]
        temp <- try( linkSpe[[k]][1] <- as.numeric(linkSpe[[k]][1]),silent = FALSE)
        if(inherits(temp ,'try-error') | temp < 2){
          stop("for I-splines link function: the first argument must be an integer greater than 1 (number of knots)")
        }
        nknots <- as.numeric(linkSpe[[k]][1])
        
        temp <- try(linkSpe[[k]][3] <- as.numeric(linkSpe[[k]][3]),silent = FALSE)
        if(inherits(temp ,'try-error') | temp < 1){
          stop("for I-splines link function: the second argument  must be a positive integer (degree of splines)")
        }
        degree[k] <- as.numeric(linkSpe[[k]][3])
        
        if(!(linkSpe[[k]][2] %in% c("quant", "manual", "equi"))){
          stop("the type of knots must be within: quant, manual or equi")
        }
        if(linkSpe[[k]][2] == "manual" & (is.null(knots[[k]]) | (length(knots[[k]])!=nknots))){
          stop("When specified manually, the number of knots must match 
               the first argument of the link function specification")
        }
        if(linkSpe[[k]][2] == "manual" & (length(knots[[k]])== nknots)){
          min <- min(Y[,col[k]], na.rm = TRUE)
          max <- max(Y[,col[k]], na.rm = TRUE)
          if( (min > min(knots[[k]])) & (max < max(knots[[k]])))stop("Knots must be in the range of the outcome")
        }
        if(linkSpe[[k]][2] == "equi" & (!is.null(knots[[k]]))){
          stop("When specified as equidistant, there is no need to specify manually the position of knots")
        }
        
        #
        if(linkSpe[[k]][2] == "quant" & (!is.null(knots[[k]]))){
          stop("When specified as placed at quantiles, there is no need to specify manually the position of knots")
        }
        if(linkSpe[[k]][2] == "quant" & (is.null(knots[[k]]))){
          knots[[k]] <- as.vector(quantile(Y[,col[k]], probs = seq(from = 0, to = 1, by = 1/(nknots-1)), na.rm = TRUE))
        }
        if(linkSpe[[k]][2] == "equi"){
          knots[[k]] <- seq(from = minY[k], to = maxY[k], by = (maxY[k]-minY[k])/(nknots-1))
        }
        
        ## if two quantiles are equal
        for(nk in 3:nknots){
          if(knots[[k]][nk]== knots[[k]][nk-1]) knots[[k]][nk-1] <- knots[[k]][nk-1] - (max(Y[,col[k]], na.rm = TRUE)-min(Y[,col[k]], na.rm = TRUE))/5
        }
        if(nknots>2){
          int_knots <- as.vector(as.numeric(knots[[k]][-c(1,nknots)]))
        }else{
          int_knots <- NULL
        }
        
        
        modISpline <- paste("~ 1 + splines2::iSpline(",col[k],",knots=","int_knots",",","degree=", degree[k],
                            ",", "intercept = T,", "derivs= 0,", "Boundary.knots= c(",minY[k],",",maxY[k],"))")
        
        modMSpline <- paste("~ -1 + splines2::iSpline(",col[k],",knots=","int_knots",",","degree=", degree[k],
                            ",", "intercept = T, ", "derivs = 1,", "Boundary.knots= c(",minY[k],",",maxY[k],"))")
        
        
        IsMat <- model.matrix(as.formula(modISpline), data = Y, na.action = na.action)
        MsMat <- model.matrix(as.formula(modMSpline), data = Y, na.action = na.action)
        colnamesY <- c(colnamesY, paste(outcomes[k],link[k], seq(1,ncol(IsMat)), sep = "."))
        colnamesYPrim <- c(colnamesYPrim, paste(outcomes[k],link[k], seq(1,ncol(MsMat)), sep = "."))
        
        Mod.MatrixY <- cbind(Mod.MatrixY, IsMat)
        Mod.MatrixYprim <- cbind(Mod.MatrixYprim, MsMat)
        df <-c(df, ncol(IsMat))
        }
  }
    Mod.MatrixY <- as.matrix(Mod.MatrixY)
    Mod.MatrixYprim <- as.matrix(Mod.MatrixYprim)
    colnames(Mod.MatrixY) <- colnamesY
    colnames(Mod.MatrixYprim) <- colnamesYPrim 
    return(list(minY = minY, maxY =maxY, degree = degree, knots=knots, df=df, Mod.MatrixY = Mod.MatrixY, 
                Mod.MatrixYprim = Mod.MatrixYprim))
    na.action = current.na.action
  }else{
    stop("Need package splines2 to work, please install it.")
  }
}

#=====================================================================================
# create object of type list containing:
# - design matrix for all sub-models of the  principal model
#====================================================================================
DataFormat <- function(data, subject, fixed_X0.models , randoms_X0.models , fixed_DeltaX.models, 
                       randoms_DeltaX.models, mod_trans.model, link = NULL, knots = NULL, 
                       outcomes, nD, Time, DeltaT){
  
  cl <- match.call()
  colnames<-colnames(data)
  id_and_Time <- data[,c(subject,Time)]
  Ni<-unique(data[,subject])
  I <- length(Ni)# number of visit
  K <- length(outcomes)
  
  # Pre-traitement of data : delete lignes with no observation
  d <- as.data.frame(data[,outcomes])
  R <- as.numeric(apply(X = d, MARGIN = 1, FUN = is_na_vec))
  data <-data[which(R==0),]
  
  
  #all predictor of the model==============================================================
  colnames <- colnames(data)
  all.pred.fixed_X0 <- NULL
  all.pred.fixed_DeltaX <- NULL
  all.pred.randoms_X0 <- NULL
  all.pred.randoms_DeltaX <- NULL
  all.pred.mod_trans <- NULL
  for( n in 1: nD){
    all.pred.fixed_X0 <- c(all.pred.fixed_X0,list(strsplit(fixed_X0.models[n],"[+]")[[1]]))
    all.pred.fixed_DeltaX <- c(all.pred.fixed_DeltaX,list(strsplit(fixed_DeltaX.models[n],"[+]")[[1]]))
    all.pred.randoms_X0 <- c(all.pred.randoms_X0,list(strsplit(randoms_X0.models[n],"[+]")[[1]]))
    all.pred.randoms_DeltaX <- c(all.pred.randoms_DeltaX,list(strsplit(randoms_DeltaX.models[n],"[+]")[[1]]))
  }
  #
  all.pred.fixed_X0<-sapply(all.pred.fixed_X0,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.fixed_X0<-sapply(all.pred.fixed_X0,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #
  all.pred.fixed_DeltaX<-sapply(all.pred.fixed_DeltaX,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.fixed_DeltaX<-sapply(all.pred.fixed_DeltaX,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #
  all.pred.randoms_X0<-sapply(all.pred.randoms_X0,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.randoms_X0<-sapply(all.pred.randoms_X0,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #
  all.pred.randoms_DeltaX<-sapply(all.pred.randoms_DeltaX,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.randoms_DeltaX<-sapply(all.pred.randoms_DeltaX,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #
  all.pred.mod_trans <- c(all.pred.mod_trans,list(strsplit(mod_trans.model,"[+]")[[1]]))
  all.pred.mod_trans<-sapply(all.pred.mod_trans,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
  all.pred.mod_trans<-sapply(all.pred.mod_trans,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
  #ajout
  all.preds<-unlist(unique(c(unlist(all.pred.fixed_X0), unlist(all.pred.fixed_DeltaX), 
                             unlist(all.pred.randoms_X0), unlist(all.pred.randoms_DeltaX),
                             all.pred.mod_trans)))
  all.preds<-inclu_intercerpt(all.preds) # remplace 1 par "(Intercept)"
  all.preds<-unique(unlist(sapply(all.preds,FUN=function(x)strsplit(x,":")[[1]])))
  
  ###all.pred_san_inter signifie all.pred_sans_intercept
  all.preds<-all.preds[-which(all.preds %in% c("(Intercept)"))]
  all.preds <- c(all.preds,Time)
  all.preds <- unique(all.preds[which(all.preds %in% colnames)])  
  
  #Case of  unobserved components  at time t
  m_i <-as.data.frame(table(as.factor(data[,subject])))$Freq # matrice of frequencies m_i
  tau_is <- data[,Time]/DeltaT # vector of individuals visits vectors
  Tmax <- max(tau_is,na.rm = TRUE)
  tau <- 0:Tmax  
  Y <- NULL
  IND <- NULL
  indY <- NULL
  data0 <- NULL
  ###cration de data0==========
  ## for x and z
  all.Y<-seq(1,K)
  for (k in 1:K)
  {
    dtemp <- data[,c(subject,outcomes[k],all.preds)]
    Y <- c(Y, dtemp[,outcomes[k]])
    IND <- c(IND, dtemp[,subject])
    indY <- c(indY,rep(all.Y[k],nrow(dtemp)))
    data0 <- rbind(data0, dtemp[,c(setdiff(colnames(dtemp),outcomes[k]))])
  }
  
  data0<-cbind(data0,Y,indY)
  data0<-data0[order(data0[,subject]),]
  data0<-na.omit(data0)
  ## for x and z
  indY <- data0[,"indY"]
  Y<-data0[,"Y"]
  
  
  #### only for x0 and x ###############
  xs <- as.data.frame(unique(data0[,c(subject,setdiff(all.preds,Time))]))
  colnames(xs) <- c(subject,setdiff(all.preds,Time))
  qsz <- as.data.frame(xs[rep(row.names(xs), rep(length(tau), dim(xs)[1])),])
  colnames(qsz) <- c(subject,setdiff(all.preds,Time))
  Times <- as.data.frame(DeltaT*rep(tau, I))
  colnames(Times) <- Time
  data_c0 <- cbind(qsz,Times)
  data_xzMatA_cov <-data_c0  
  data_xzMatA_cov <-data_xzMatA_cov[order(data_xzMatA_cov[,subject], data_xzMatA_cov[,Time]),]
  data_xzMatA_cov <- data_xzMatA_cov[!duplicated(data_xzMatA_cov),]
  x_cov<- NULL
  
  for(n in 1:nD){
    indLP_x <- rep(n, dim(data_xzMatA_cov)[1])
    data_x_cov_i <- cbind(data_xzMatA_cov, indLP_x)
    x_cov <- rbind(x_cov, data_x_cov_i)
  }
  x_cov <- x_cov[order(x_cov[,subject],x_cov[,Time]),]
  
  
  ##only for x0 #####
  x0 <- NULL
  nb_x0_n <- NULL
  col_n<-list()
  x0_cov <- x_cov[which(x_cov[,Time]==0),]
  #
  ##
  indLP_x0 <- x0_cov$indLP_x
  for(n in 1:nD){
    r<-as.formula(paste(subject, fixed_X0.models[n], sep="~-1+"))
    x0n<-model.matrix(r,data=x0_cov)
    nb_x0_n <- c(nb_x0_n,ncol(x0n))
    if(length(x0n)==0){
      col <- paste(n,"zero",sep="")
      x0n<-matrix(assign(col,rep(0,dim(x0_cov)[1])))
      nb_x0_n <- c(nb_x0_n,ncol(x0n))
    }
    
    colnames<-colnames(x0n)
    colnames<-paste("LP0",n,colnames,sep=".")
    colnames(x0n) <-colnames
    col_n <-c(col_n,list(colnames))
    x0<-cbind(x0,as.matrix(x0n))
  }
  x0 <-cbind(indLP_x0,x0)
  ### remplissage avec les zeros
  tous_col_x0 <-unlist(col_n)
  for(i in 1:nrow(x0)){
    col_i <- unlist(col_n[[x0[i,"indLP_x0"]]])
    col_0<-tous_col_x0[which(!(tous_col_x0 %in% col_i))]
    x0[i,col_0]<-0 #  passer pour optimisation
  }
  x0 <- as.matrix(x0)
  colnames <- colnames(x0)
  x0 <- as.matrix(x0[,-c(1)])
  colnames(x0) <- colnames[-c(1)]
  #   x0 <- as.matrix(x0)
  
  
  ##only for x #####
  x <- NULL
  nb_x_n <- NULL
  col_n<-list()
  indLP_x <- x_cov$indLP_x
  for(n in 1:nD){
    r<-as.formula(paste(subject, fixed_DeltaX.models[n], sep="~-1+"))
    xn<-model.matrix(r,data=x_cov)
    nb_x_n <- c(nb_x_n,ncol(xn))
    if(length(xn)==0){
      col <- paste(n,"zero",sep="")
      xn<-matrix(assign(col,rep(0,dim(x_cov)[1])))
      nb_x_n <- c(nb_x_n,ncol(xn))
    }
    colnames<-colnames(xn)
    colnames<-paste("DeltaLP",n,colnames,sep=".")
    colnames(xn) <-colnames
    col_n <-c(col_n,list(colnames))
    x<-cbind(x,as.matrix(xn))
  }
  x <-cbind(indLP_x,x)
  ### filling with zeros
  tous_col_x <-unlist(col_n)
  for(i in 1:nrow(x)){
    col_i <- unlist(col_n[[x[i,"indLP_x"]]])
    col_0<-tous_col_x[which(!(tous_col_x %in% col_i))]
    x[i,col_0]<-0 # z  passer pour optimisation
  }
  
  x <- as.matrix(x)
  colnames <- colnames(x)
  x <- as.matrix(x[,-c(1)])
  colnames(x) <- colnames[-c(1)]
  
  #===================================================================
  #     construction  of matrices z0 et z========================
  data_z_cov <- data_xzMatA_cov[, c(subject,Time)]
  z_cov<- NULL
  for(n in 1:nD){
    indY_z <- rep(n, dim(data_z_cov)[1])
    data_z_cov_i <- cbind(data_z_cov, indY_z)
    z_cov <- rbind(z_cov, data_z_cov_i)
  }
  z_cov <- z_cov[order(z_cov[,subject],z_cov[,Time]),]
  #   z_cov[,Time] <- z_cov[,Time]*DeltaT ######
  
  #### only for z0 ####
  z0_cov <- z_cov[which(z_cov[,Time]==0),]
  indY_z0 <- z0_cov$indY_z
  z0 <- NULL
  col_n<-list()
  q0 <- NULL
  nb_paraDw <- 0
  for(n in 1:nD){
    r<-as.formula(paste(subject,randoms_X0.models[n], sep="~-1+"))
    z0n<-model.matrix(r,data=z0_cov)
    if(length(z0n)==0){
      col <- paste(n,"zero",sep="")
      z0n<-matrix(assign(col,rep(0,dim(z0_cov)[1])))
    }
    colnames<-colnames(z0n)
    colnames<-paste(n,colnames,sep="")
    colnames(z0n) <-colnames
    col_n <-c(col_n,list(colnames))
    z0<-cbind(z0,z0n)
    q0 <- c(q0,ncol(z0n))
    
    
  }
  z0 <-cbind(indY_z0,z0)
  ### filling with zeros
  tous_col_z0 <-unlist(col_n)
  for(i in 1:nrow(z0)){
    col_i <- unlist(col_n[[z0[i,"indY_z0"]]])
    col_0<-tous_col_z0[which(!(tous_col_z0 %in% col_i))]
    z0[i,col_0]<-0 # z  passer pour optimisation
  }
  z0 <- z0[,-c(1)]
  z0 <- as.matrix(z0)
  
  
  #### only for z ####
  indY_z <- z_cov$indY_z
  z <- NULL
  col_n<-list()
  q <- NULL
  nb_paraDu <- 0
  for(n in 1:nD){
    r<-as.formula(paste(subject,randoms_DeltaX.models[n], sep="~-1+"))
    zn<-model.matrix(r,data=z_cov)
    if(length(zn)==0){
      col <- paste(n,"zero",sep="")
      zn<-matrix(assign(col,rep(0,dim(z_cov)[1])))
    }
    colnames<-colnames(zn)
    colnames<-paste(n,colnames,sep="")
    colnames(zn) <-colnames
    col_n <-c(col_n,list(colnames))
    z<-cbind(z,zn)
    q <- c(q,ncol(zn))
    
  }
  z <-cbind(indY_z,z)
  ### filling with zeros
  tous_col_z <-unlist(col_n)
  for(i in 1:nrow(z)){
    col_i <- unlist(col_n[[z[i,"indY_z"]]])
    col_0<-tous_col_z[which(!(tous_col_z %in% col_i))]
    z[i,col_0]<-0 #
  }
  z <- z[,-c(1)]
  z <- as.matrix(z)
  #============================================================
  # design matrix for transition model
  f<-as.formula(paste(subject,mod_trans.model, sep="~"))# put subject, just to have a left side for the formula
  
  modA_mat<-model.matrix(as.formula(paste(subject,mod_trans.model, sep="~")),data=data_xzMatA_cov)
  #   dim(modA_mat)
  #   head(modA_mat)
  #============================================================
  #design matrix for markers transformation
  Y <- as.matrix(data[,outcomes])
  tr_Y <- f.link(outcomes = outcomes, Y=as.data.frame(Y), link=link, knots =knots)
  Mod.MatrixY <- tr_Y$Mod.MatrixY
  Mod.MatrixYprim <- tr_Y$Mod.MatrixYprim
  knots <- tr_Y$knots
  degree = tr_Y$degree
  df <- tr_Y$df
  minY <- tr_Y$minY
  maxY <- tr_Y$maxY
  nb_RE <- sum(q0,q)
  nb_paraD <- nb_RE*(nb_RE+1)/2
  
  return(list(nb_subject=I, nb_obs = length(na.omit(as.vector(Y))), K=K, nD = nD, all.preds = all.preds, id_and_Time=id_and_Time,Tmax = Tmax, m_i = m_i, Y = Y, Mod.MatrixY=Mod.MatrixY,  
              Mod.MatrixYprim=Mod.MatrixYprim, minY = minY, maxY = maxY, knots = knots, df = df, degree = degree, x = x, x0 = x0, 
              vec_ncol_x0n = nb_x0_n, z = z, z0=z0, q = q, q0 = q0, nb_paraD = nb_paraD, nb_RE=nb_RE, modA_mat = modA_mat,
              tau = tau, tau_is = tau_is))
}







