% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{tsGmatprodA0totau}
\alias{tsGmatprodA0totau}
\title{Function that  creates a big matrix ts_G_mat_prod_A_0_to_tau containing  Prod(A_t)t=0,tau.}
\usage{
tsGmatprodA0totau(K, vec_alpha_ij, tau, DeltaT, modA_mat)
}
\arguments{
\item{K}{an integer representing the size of K*K matrix}

\item{vec_alpha_ij}{a vector of overall parameters associated to the
model.matrix for elements of the transistion matrix}

\item{tau}{vector of integers indicating times}

\item{DeltaT}{double that indicates the discretization step}

\item{modA_mat}{model.matrix for elements of the transistion matrix}
}
\value{
a matrix
}
\description{
Function that  creates a big matrix ts_G_mat_prod_A_0_to_tau containing  Prod(A_t)t=0,tau.
}
