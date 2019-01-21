#==================================================================================================
#' Function to plot trajectories of observation and their predictions
#'
#' @param Obs Observations
#' @param outcomes names of outcomes
#' @param preds indicates predictions
#' @param VarTime indicates the name of the covariate representing the time
#' @param with_IC indicates if the plot is done with confidence band, default value is FALSE
#' @param Ylim minimum and maxmum of Y axe
#' @param xlab xlab
#' @param ylab ylab
#' @param age0 indicates baseline age
#' @param Title Titles for the plots
#'
#' @return ---
#' @export


PlotX <- function(Obs, outcomes, preds, VarTime, with_IC=FALSE, Ylim,  xlab = "Time", ylab = "Score", age0=0, Title){
  # k designe le marqueur 
  K <- length(outcomes)
  preds.colnames <- colnames(preds)
  Pred_mean <- NULL
  Time <- sort(unique(preds[,VarTime]))
  for(t in 1: length(Time)){
    Pred_mean_t <- apply(preds[which(preds[,VarTime]==Time[t]),],MARGIN = 2, FUN = mean, na.rm =TRUE)
    Pred_mean <- rbind(Pred_mean,Pred_mean_t)
  }
  Pred_mean <- data.frame(Pred_mean)
  # Pred_mean$t <- age0 + Pred_mean$t
  # par(mfrow = c(round(K/2),round(K/2)))
  #     par(mfrow = c(round(K/2),round(K/2)), cex = 1, mar=c(5,5,2,2) )
  for(k in 1:K){
    b_inf <- Obs$obs[which(!is.nan(Obs$obs[,outcomes[k]])),outcomes[k]]-1.96*Obs$se_obs[which(!is.nan(Obs$obs[,outcomes[k]])),outcomes[k]]
    b_sup <- Obs$obs[which(!is.nan(Obs$obs[,outcomes[k]])),outcomes[k]]+1.96*Obs$se_obs[which(!is.nan(Obs$obs[,outcomes[k]])),outcomes[k]]
    
    # Ylim = c(min(Obs$obs[,outcomes[k]],b_inf, na.rm = TRUE), max(Obs$obs[,outcomes[k]],b_sup, na.rm = TRUE))
    
    plot(y=Obs$obs[which(!is.nan(Obs$obs[,outcomes[k]])),outcomes[k]], x=(Obs$obs[which(!is.nan(Obs$obs[,outcomes[k]])),1]),
         type = "l",pch = 20, col =1, ylim = Ylim[[k]], lwd=3, xaxt = "n", xlab = xlab, ylab = ylab, cex.lab = 1.5
    )
    #     axis(1, (age0+Obs$obs[,1]*pas), (age0+Obs$obs[,1]*pas),1)
    #     axis(2)
    axis(1, seq(from = 0, to = 4, by = 0.23), seq(from = 0, to = 4, by = 0.23),1)
    axis(2)
    
    if(with_IC==TRUE){
      lines(y=b_inf, x=(Obs$obs[which(!is.nan(Obs$obs[,outcomes[k]])),1]), type = "b", pch = -1, col =1, lty = 2, lwd=2.5)
      lines(y=b_sup, x=(Obs$obs[which(!is.nan(Obs$obs[,outcomes[k]])),1]), type = "b", pch = -1, col =1, lty = 2, lwd=2.5)
    }
    #Prediction
    #lines(Pred_mean[which(!is.nan(Obs$obs[,(k+1)])),(k+1)]~Pred_mean$t[which(!is.nan(Obs$obs[,(k+1)]))], type = "l", col ="grey", lty=2, lwd=4) # prediction marginale
    #lines(Pred_mean[which(!is.nan(Obs$obs[,(k+1)])),(k+K+1)]~Pred_mean$t[which(!is.nan(Obs$obs[,(k+1)]))], type = "b", lty=2, col ="grey", lwd=4)
    #     points(Pred_mean[,c(1,(k+1))], pch = 17, col ="grey", lwd=3)
    points(Pred_mean[,c(which(preds.colnames==VarTime),which(preds.colnames == paste(outcomes[k],sep = "-","Pred")))], type = "p", pch = 4, col ="black", cex = 1, lwd = 3)
    #     legend <- c("observations", "band confidence of obs.", "marginal predictions", "subject specific predictions")
    #     legend(x=pos.legend, legend,  bty = "n", cex = .5, lwd = 0.5, col=c("black","black", "grey","grey"), lty=c(2,2,1,2), pch=c(4,-1,-1,16)) # puts text in the legend
    title(main = Title[k], cex.main = 1.5)
  }
}
