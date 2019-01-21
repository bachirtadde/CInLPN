#' Generate a multivate normal vector
#'
#' @param seed seed for the randomization
#' @param m mean of distribution
#' @param sd standard deviation of the distribution
#'
#' @return a multivate normal vector

f_mvrnorm<- function(seed, m, sd){
  if(requireNamespace("MASS", quietly = TRUE)){
    set.seed(seed)
    return(MASS::mvrnorm(n=1, mu = m, Sigma = sd))
  }else{
    stop("Need package MASS to work, Please install it.")
  }
}
