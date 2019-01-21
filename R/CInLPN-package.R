#' A package for Causal Inference in a Latent Processes Network
#'
#' This package Contains functions to fit a dynamic model for multiple latent processes. 
#' The methodology is based on the combinaison 
#' of a linear mixed model for the latent processes trajectories and a system of 
#' differences equations for assessing their temporal influences. 
#' The estimation is done in the maximum likelihood framework.
#'
#' \tabular{ll}{ Package: \tab CInLPN\cr Type: \tab Package\cr
#' Version: \tab 0.2.0\cr Date: \tab 2019-01-03\cr License: \tab GPL (>= 2.0)\cr}
#' 
#' @exportPattern ^[[:alpha:]]+
#' @importFrom graphics abline axis lines par plot points title
#' @importFrom stats as.formula model.matrix na.action na.omit pchisq pnorm printCoefmat quantile terms var
#' @importFrom Rcpp evalCpp 
#' @name CInLPN-package
#' @docType package
#' @author Bachirou Tadd\'e, C\'ecile Proust-Lima
#' Tadd\'e O. B., et al. (2018). Dynamic modeling of multivariate latent processes 
#' and their temporal relationships:  Application to Alzheimer's disease,  
#' (submitted paper available at https://arxiv.org/abs/1806.03659)
#' 
#' @keywords "Causality", "Dynamic model"," Latent processes"," multivariate longitudinal data"
#' @useDynLib CInLPN, .registration = TRUE
NULL
