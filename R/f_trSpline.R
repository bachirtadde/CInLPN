#' Function that compute  transformaitions  by I-splines basis and their derivates
#'
#' @param y outcomes in real scales
#' @param minY minimum of outcome y
#' @param maxY maximum of outcome y
#' @param knots indicates position of knots used to transform outcomes
#' @param degree indicates degree of the basis of splines
#' @param paras initial values for parameters
#'
#' @return a matrix

f_trSpline <- function(y, minY, maxY, knots, degree, paras){
  if(requireNamespace("splines2", quietly = TRUE)){
    ni <- length(y)
    y <- data.frame(y=y)
    nk <- length(knots)
    int_knots <- knots[-c(1,nk)]
    
    modISpline <- paste("~ 1 + splines2::iSpline( y", ",knots=","int_knots",",","degree=", degree,
                        ",", "intercept = T,", "derivs= 0,", "Boundary.knots= c(",minY,",",maxY,"))")
    modMSpline <- paste("~ -1 + splines2::iSpline(y",",knots=","int_knots",",","degree=", degree,
                        ",", "intercept = T, ", "derivs= 1,", "Boundary.knots= c(",minY,",",maxY,"))")
    
    IsMat <- model.matrix(as.formula(modISpline), data = y, na.action = na.action)
    MsMat <- model.matrix(as.formula(modMSpline), data = y, na.action = na.action)
    
    colnames(IsMat) <- paste("y", seq(1,ncol(IsMat)), sep = ".")
    colnames(MsMat) <- paste("y", seq(1,ncol(MsMat)), sep = ".")
    tr <- matrix(rep(NA,ni*2), ncol = 2, byrow = F)
    tr[,1] <- IsMat%*%paras
    tr[,2] <- 1/(MsMat%*%paras[-1])
    colnames(tr) <- c("Y.tr","invJ.Y.tr")
    return(tr)
  }else{
    stop("Need package splines2 to work, please install it.")
  }
}
