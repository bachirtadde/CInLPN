f_mvrnorm<- function(seed, m, sd){
  if(requireNamespace("MASS", quietly = TRUE)){
    set.seed(seed)
    return(MASS::mvrnorm(n=1, mu = m, Sigma = sd))
  }else{
    stop("Need package MASS to work, Please install it.")
  }
}
