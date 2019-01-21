#' Summary of CInLPN object
#'
#' @param object an CInLPN object
#' @param \dots optional parameters
#'
#' @return 0
#' @export
summary.CInLPN <- function(object, ...){
  x <- object
  if (!inherits(x, "CInLPN")) stop("use only with \"CInLPN\" objects")
  
  # cat("Dynamic Temporal links between latent processes", "\n")
  
  
  cl <- x$call
  
  cat("Model fitted by maximum likelihood", "\n")
  
  cl$B <- NULL
  if(is.data.frame(cl$data))
  {
    cl$data <- NULL
    x$call$data <- NULL    
  }
  cat(" \n")
  dput(cl)
  cat(" \n")
  
  # posfix <- eval(cl$posfix)
  
  cat("Statistical Model:", "\n")
  cat(paste("     Dataset:", x$call$data),"\n")
  cat(paste("     Number of subjects:", x$ns),"\n")
  cat(paste("     Number of observations:", x$N),"\n")
  #if(length(x$na.action))cat(paste("     Number of observations deleted:",length(x$na.action)),"\n")
  cat(paste("     Number of latent process(es):", x$nD), "\n")
  cat(paste("     Number of marker(s):", x$K), "\n")
  cat(paste("     Number of parameters:", length(x$best))," \n")
  #if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")
  cat("     Link function(s): ")
  for (yk in 1:x$K)
  {
    if (is.null(x$linkstype[yk]) || (gsub("[[:space:]]","",x$linkstype[yk])=="linear")) {
      if (yk>1) cat("                       ")
      cat("linear for",x$outcomes[yk]," \n")
    }
    else{
      if (yk>1) cat("                       ")
      cat(gsub("[[:space:]]","",x$linkstype[yk]),  " at knodes:", round(x$linknodes[[yk]],3)," for ",x$outcomes[yk], "\n")
    }
  }
  
  cat(" \n")
  cat("Optimisation process:", "\n")
  
  if(x$conv==1) cat("     Convergence criteria satisfied")
  if(x$conv==2) cat("     CAUTION: Maximum number of iteration reached without convergence")
  if(x$conv==3) cat("     Convergence with restrained Hessian matrix")
  if(x$conv==4|x$conv==12) {
    cat("     The program stopped abnormally. No results can be displayed.\n")
  }
  else{
    
    cat(" \n")
    cat("     Number of iterations: ", x$niter, "\n")
    cat("     Convergence criteria: parameters=", signif(x$ca,2), "\n")
    cat("                         : likelihood=", signif(x$cb,2), "\n")
    cat("                         : second derivatives=", signif(x$rdm,2), "\n")
    cat(" \n")
    cat("Goodness-of-fit statistics:", "\n")
    cat(paste("     maximum log-likelihood:", round(x$loglik,2))," \n")
    cat(paste("     AIC:", round(x$AIC,2))," \n")
    cat(paste("     BIC:", round(x$BIC,2))," \n")
    cat(" \n")
    
    cat("Maximum Likelihood Estimates:", "\n")
    cat(" \n")
    
    m <- length(x$best)
    var_cov <- matrix(0,nrow=m,ncol=m)
    var_cov[upper.tri(var_cov,diag=TRUE)] <- x$v
    var_cov[lower.tri(var_cov,diag=FALSE)] <- t(var_cov)[lower.tri(var_cov,diag=FALSE)]
    
    #
    se <- rep(NA,length(x$coefficients))
    tval <- rep(NA,length(x$coefficients))
    p.value <- rep(NA,length(x$coefficients))
    se[which(x$posfix==0)] <-sqrt(diag(var_cov))
    #test
    tval[which(x$posfix==0)] <- x$b/ na.omit(se)
    p.value[which(x$posfix==0)] <- 2*pnorm(abs(na.omit(tval)),mean = 0, sd = 1, lower.tail = FALSE)
    
    TAB <- cbind(Estimate = x$coefficients,
                 StdErr = se,
                 Wald = tval,
                 p.value = p.value)
    # printCoefmat(cbind(num=seq(from = 1, to = length(x$coefficients), by = 1),round(TAB,4)), P.value=TRUE, has.Pvalue=TRUE)
    printCoefmat(round(TAB,4), P.values =TRUE, has.Pvalue=TRUE)
    
  }
}