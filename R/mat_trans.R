### function to create a matrix of all transition matrix and their standard deviation
mat_trans <- function(model, K, data, ind, DeltaT, ind_parafixA, val_parafixA){
  cl <- match.call()
  
  # if(requireNamespace("splines", quietly = TRUE)){
  ##
  nb_para_mu0 <- model$length_para_mu0
  nb_para_mu <- model$length_para_mu
  nb_para_RE <- model$length_para_RE
  vec_alpha_ij <- rep(0, length(ind_parafixA))
  vec_alpha_ij[which(ind_parafixA==0)]<- model$best[(nb_para_mu0+nb_para_mu+nb_para_RE+1):
                                                      (nb_para_mu0+nb_para_mu+nb_para_RE+length(ind_parafixA[which(ind_parafixA==0)]))]
  vec_alpha_ij[which(ind_parafixA==1)]<- val_parafixA
  nb_para_a <- 0
  if(!is.null(ind_parafixA[which((ind_parafixA==0))])){
    nb_para_a <- length(ind_parafixA[which((ind_parafixA==0))])
  }
  
  
  #
  maxTime <- (max(model$tau)-1)*DeltaT
  Time <- seq(from = 0, to = maxTime, by = 0.1)
  tau <- Time/0.1
  ntau <- length(tau)
  matA <- matrix(rep(0, ntau*K*K), ncol = K*K)
  SeA <- matrix(rep(0, ntau*K*K), ncol = K*K)
  #============================================================
  ## computing of design matrix 
  
  I <- dim(data)[1]
  data_matA <- NULL
  for(i in 1:I){
    data_mat_i <- NULL
    for(t in tau){
      data_mat_i <- rbind(data_mat_i, data[i,])
    }
    data_matA <- rbind(data_matA, data_mat_i)
  }
  
  Time <- data.frame(Time)
  colnames(Time)<- model$Time
  data_matA <- data.frame(data_matA,Time)
  
  
  f<-as.formula(model$call$structural.model$trans.matrix)# nom.subject pour juste avoir un premier membre pour la formule
  modA_mat<-model.matrix(f,data=data_matA)
  
  df <- ncol(modA_mat)
  # Compute variance-covariance matrix====
  m <- length(model$best)
  var_cov <- matrix(0,nrow=m,ncol=m)
  var_cov[upper.tri(var_cov,diag=TRUE)] <- model$v
  var_cov[lower.tri(var_cov,diag=FALSE)] <- t(var_cov)[lower.tri(var_cov,diag=FALSE)]   
  # var_cov <-res$Var
  var_cov.vec_alpha_ij <- matrix(rep(0,length(ind_parafixA)*length(ind_parafixA)), ncol = length(ind_parafixA))
  if(nb_para_a!=0){
    vA.est <- as.matrix(var_cov[(nb_para_mu0+nb_para_mu+nb_para_RE+1):(nb_para_mu0+nb_para_mu+nb_para_RE+nb_para_a),
                                (nb_para_mu0+nb_para_mu+nb_para_RE+1):(nb_para_mu0+nb_para_mu+nb_para_RE+nb_para_a)])
    vecA <- as.vector(vectorise(vA.est)) # vectorsise
  }
  iA <- 1
  for(ii in 1:length(ind_parafixA)){
    for(jj in 1:length(ind_parafixA)){
      if(ind_parafixA[ii] == 0 & ind_parafixA[jj]==0){
        var_cov.vec_alpha_ij[ii,jj] <- vecA[iA]
        var_cov.vec_alpha_ij[jj,ii] <- vecA[iA]
        iA <- iA+1
      }
    }
  }
  
  # initialisation======
  matA <- matrix(rep(0, ntau*K*K), ncol = K*K)
  SeA <- matrix(rep(0, ntau*K*K), ncol = K*K)
  modA_mat_i <- modA_mat[((ind-1)*ntau+1):((ind)*ntau),]
  for(t in 1: ntau){
    matA[t,] <- vecaijt( K, (t-1), vec_alpha_ij, as.matrix(modA_mat_i))
  }
  #compute de Se_a, the standard deviation
  p=0
  for(k in 1:(K*K)){
    var_cov.vec_alpha_a_i <- as.matrix(var_cov.vec_alpha_ij[(p+1):(p+df),(p+1):(p+df)])
    SeA[,k] <- sqrt(abs(diag(modA_mat[((ind-1)*ntau+1):(ind*ntau),]%*%var_cov.vec_alpha_a_i%*%t(modA_mat[((ind-1)*ntau+1):(ind*ntau),]))))
    p <- p+df
  } 
  
  titre_a <-NULL
  titre_a1 <- seq(1,K)
  titre_se <-NULL
  titre_se1 <- seq(1,K)
  for(k in 1:K){
    titre_a<- c(titre_a, paste(titre_a1[k], titre_a1,  sep=""))
    titre_se<- c(titre_se, paste(titre_se1[k], titre_se1,  sep=""))
  }
  titre_a <- paste("a_", titre_a, sep="")
  titre_se <- paste("se_", titre_se, sep="")
  
  colnames(matA) <- titre_a
  matA <- as.data.frame(cbind(Time, matA))
  colnames(SeA) <- titre_se
  SeA <- as.data.frame(cbind(Time, SeA))
  return(list(matA=matA, SeA=SeA))
  # }else{
  #   stop("Need package splines to fit model using splines basis to model temporal influences. \n Please install
  #        splines package before.")
  # }
}
