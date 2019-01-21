#' Function used to replace 1 by intercept in a formula
#'
#' @param x vector of text
#'
#' @return vector of characters

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
#' Function that compute the length of vector without NAs
#'
#' @param x numeric vector with possibly NAs
#'
#' @return a numeric vector

taille <- function(x){
  x <- na.omit(x)
  return(length(x))
}




#' function to create the variance-covariance matrix  from vector of all element of the matrix
#'
#' @param model.b estimates of parameters
#' @param model.v variance-covariance elements
#'
#' @return a matrix

f_var_cov <- function(model.b,model.v){
  cl <- match.call()
  m <- length(model.b)
  var_cov <- matrix(0,nrow=m,ncol=m)
  var_cov[upper.tri(var_cov,diag=TRUE)] <- model.v
  var_cov[lower.tri(var_cov,diag=FALSE)] <- t(var_cov)[lower.tri(var_cov,diag=FALSE)]
  var_cov
}