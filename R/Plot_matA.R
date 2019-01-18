#=========== plot
Plot_a <- function(recaps.A, with_IC=FALSE, Ylim=c(-2,2), age0, legend, pos.legend, color, xlab ="Time (in year)"){ 
  pas = 0.1
  Nb.model <- length(recaps.A) ## nombre de modeles
  if(length(age0)==1) age0 <- rep(age0,Nb.model)
  if(length(with_IC)==1) with_IC <- rep(with_IC,Nb.model)
  if(Nb.model < 1) stop("at least one model is needed")
  ## check for convergence of all models
  for(n in 1: Nb.model){
    if(is.null(recaps.A[[n]]))stop("one model doesn't reach convergence: impossible to plot matrix A")
  }
  K <- sqrt(dim(recaps.A[[1]]$matA)[2]) ## nombre de spheres
  
  ## cols des graphiques
  cols <-colnames(recaps.A[[1]]$matA)
  # par(mfrow = c(1,1), lwd=1)
  ##====
  k <- 0
  kp <- 1
  par(mfrow = c(K,K), lwd = 2.5, cex = .5, mar=c(4,4,4,4), mai = c(0.5,0.5,0.5,0.5) )
  for(ki in 1:K){
    for(kj in 1:K){
      k <- k+1
      if(is.null(Ylim)){
        Ylim <- c(100000,-100000)
        for(n in 1: Nb.model){
          Ylim[1] <- min(Ylim[1],(min(recaps.A[[n]]$matA[,k+1])))
          Ylim[2] <- max(Ylim[2],(max(recaps.A[[n]]$matA[,k+1])))
          if( with_IC[n]==TRUE){
            Ylim[1] <- min(Ylim[1],(min(recaps.A[[n]]$matA[,k+1]-1.96*recaps.A[[n]]$SeA[,k+1])))
            Ylim[2] <- max(Ylim[2],(max(recaps.A[[n]]$matA[,k+1]+1.96*recaps.A[[n]]$SeA[,k+1])))
          }
        }
        # Ylim[1] <- max(ylim[1],Ylim[1])
        # Ylim[2] <- min(ylim[2],Ylim[2])
      }
      for(n in 1: Nb.model){
        if(n ==1){
          plot(x = age0[n]+recaps.A[[n]]$matA[,1], y =recaps.A[[n]]$matA[,k+1], type = 'l', col= color[n],lty = 1 ,ylim = Ylim,
               main = cols[k+1], xlab = xlab, ylab = "intensity", cex.main = 1.8, cex.lab=1.8, cex.axis = 1.2)
          #           Tm <- (dim(recaps.A[[1]]$matA)[1]-1)/10
          # #           axis(1, at =seq(age0[n], (age0[n]+Tm),pas))
          #           axis(1, labels =seq(age0[n], (age0[n]+Tm),pas))
          #           axis(2)
          
          if( with_IC[n]==TRUE){
            lines(x = age0[n]+recaps.A[[n]]$matA[,1], y= (recaps.A[[n]]$matA[,k+1]-1.96*recaps.A[[n]]$SeA[,k+1]), type = "l", col= color[n], lty = 2, cex.lab=1.5)
            lines(x = age0[n]+recaps.A[[n]]$matA[,1], y= (recaps.A[[n]]$matA[,k+1]+1.96*recaps.A[[n]]$SeA[,k+1]), type = "l", col= color[n], lty = 2, cex.lab=1.5 )
          }
          abline(h=0, lty=1, lwd=0.5, col = "black")
        }
        
        if(n != 1){
          lines(x = age0[n]+recaps.A[[n]]$matA[,1], y =recaps.A[[n]]$matA[,k+1], type = 'l', col= color[n],lty = 1, cex.lab=1.5)
          if( with_IC[n]==TRUE){
            lines(x = age0[n]+recaps.A[[n]]$matA[,1], y= (recaps.A[[n]]$matA[,k+1]-1.96*recaps.A[[n]]$SeA[,k+1]), type = "l", col=color[n], lty = 2)
            lines(x = age0[n]+recaps.A[[n]]$matA[,1], y= (recaps.A[[n]]$matA[,k+1]+1.96*recaps.A[[n]]$SeA[,k+1]), type = "l", col=color[n], lty = 2)
          }
          abline(h=0, lty=1, lwd=0.5, col = "black")
        }
      }
      legend(x=pos.legend[k], legend, cex = 1.5, lwd = 2.5, col=color, bty = "n", xjust = 0) # puts text in the legend
    }
    kp <- kp+K+1
  }
}