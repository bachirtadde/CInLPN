WaldMult <- function( v,beta,pos,matR, value){
  # Compute variance-covariance matrix of estimates parameters====
  # v = vector of all unique element to create the symetric matric
  m <- length(beta)
  var_cov <- matrix(0,nrow=m,ncol=m)
  var_cov[upper.tri(var_cov,diag=TRUE)] <- v
  var_cov[lower.tri(var_cov,diag=FALSE)] <- t(var_cov)[lower.tri(var_cov,diag=FALSE)]
  
 
  Mat <- matrix(0,nrow=length(pos),ncol=length(pos))
  
  Mat <- matR%*%var_cov[pos,pos]%*%t(matR)
  

  Vect <- matR%*%beta[pos]-value
  
  Wald <- t(Vect)%*%solve(Mat)%*%Vect
  
  # degre of freedom
  ddl <- length(pos)
  # Pvalue
  p_value <- 1-pchisq(Wald,df=ddl)
  return(p_value)
}
