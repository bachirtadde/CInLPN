#================================================================================================
#Function used to replace 1 by intercept in a formula
inclu_intercerpt <- function(x){
  N <- length(x)
  for(i in 1: N){
    if(x[i]=="1") x[i] <- "(Intercept)"
    if(x[i]=="-1") x[i] <- "-(Intercept)"
  }
  if("(Intercept)"%in% x & "-(Intercept)" %in% x) x <- x[-which(x=="-(Intercept)")]
  return(x)
}


## function===============
taille <- function(x){
  x <- na.omit(x)
  return(length(x))
}

# function to compute empirical expectation and empirical variance
Expect_var <- function(data, VarTime, outcomes, age0=0, pas =1){
  cl <- match.call()  
  obs <-NULL
  se_obs <- NULL
  Time <- sort(unique(data[,VarTime]))
  for(t in 1: length(Time)){
    obs_t <- apply(data[which(data[,VarTime]==Time[t]),outcomes],
                   MARGIN = 2, FUN = mean, na.rm =TRUE
    )
    obs <- rbind(obs,obs_t)
    #var
    v <- apply(data[which(data[,VarTime]==Time[t]), outcomes],MARGIN = 2, 
               FUN = var, na.rm =TRUE)
    n <- apply(data[which(data[,VarTime]==Time[t]), outcomes],MARGIN = 2, 
               FUN = taille)
    se_obs_t <- sqrt(v/n)
    se_obs <- rbind(se_obs,se_obs_t)
  }
  obs <- as.data.frame(cbind(t=Time,obs),
                       row.names =TRUE)
  se_obs <- as.data.frame(cbind(Time,se_obs),
                          row.names =TRUE)
  return(list(obs=obs, se_obs=se_obs))
}
## function to plot observed markers with confidence band at 95% ==============================
Plot_obs_Markers <- function(Obs,age0=0, pas=1, outcomes){
  cl <- match.call()
  K <- length(outcomes)
  par(mfrow = c(1,1))
  for(k in 1:K){
    #==========
    b_inf <- Obs$obs[which(!is.nan(Obs$obs[,(k+1)])),(k+1)]-1.96*Obs$se_obs[which(!is.nan(Obs$obs[,(k+1)])),(k+1)]
    b_sup <- Obs$obs[which(!is.nan(Obs$obs[,(k+1)])),(k+1)]+1.96*Obs$se_obs[which(!is.nan(Obs$obs[,(k+1)])),(k+1)]
    Ylim = c(min(Obs$obs[,(k+1)],b_inf, na.rm = TRUE), max(Obs$obs[,(k+1)],b_sup, na.rm = TRUE))
    #     plot(Obs$obs[,(k+1)]~(age0+Obs$obs[,1]*pas), type = "l",pch = 16, col =1, ylim = Ylim, 
    #          xlab = "delai", ylab = outcomes[k],cex.lab=1.5)
    plot(y=Obs$obs[which(!is.nan(Obs$obs[,(k+1)])),(k+1)], x=(age0+Obs$obs[which(!is.nan(Obs$obs[,(k+1)])),1]*pas), type = "b",pch = 16, col =1, ylim = Ylim, 
         xlab = "Time", ylab = outcomes[k],cex.lab=1.5
    )
    axis(1, (age0+Obs$obs[,1]*pas), (age0+Obs$obs[,1]*pas),1)
    axis(2)
    lines(y=b_inf, x=(age0+Obs$obs[which(!is.nan(Obs$obs[,(k+1)])),1]*pas), type = "l", col =2, lty = 2)
    lines(y=b_sup, x=(age0+Obs$obs[which(!is.nan(Obs$obs[,(k+1)])),1]*pas), type = "l", col =2, lty = 2)
  }
}

##------------------------
## function to create the variance-covariance matrix  from vector of all element of the matrix
f_var_cov <- function(model.b,model.v){
  cl <- match.call()
  m <- length(model.b)
  var_cov <- matrix(0,nrow=m,ncol=m)
  var_cov[upper.tri(var_cov,diag=TRUE)] <- model.v
  var_cov[lower.tri(var_cov,diag=FALSE)] <- t(var_cov)[lower.tri(var_cov,diag=FALSE)]
  var_cov
}

doubleAreEqual <- function(a,b){
  res <- FALSE
  if(length(a)!=length(b)) res <- FALSE
  else{
    r <- rep(NA,length(a))
    for(i in 1:length(a)){
      r[i] <- areEqual(a[i],b[i])
    }
    res <- all(r)
  }
  return(res)
}