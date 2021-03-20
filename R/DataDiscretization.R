#' function that return 1 for a vector without observation after/during discretization
#'
#' @param x a vector
#'
#' @return a binary

AnyObsByLine <- function(x){
  d<-0 
  len <-length(na.omit(x))
  if(len == 0){
    d <- 1
  }
  return(d)
}

#=====================================================================================
#' function that filters initial data to return only individuals having at least one observation for outcomes onf interest.
#' Meaning that the returned data contained individuals with at least one measurement at the considered visit time.
#'
#' @param data input data
#' @param subject subject identifiant
#' @param outcomes a vector of outcomes of interest: markers
#'
#' @return a dataframe
#' 
OneOutcomeAtLeast <- function(data, subject, outcomes){
  cl <- match.call()
  colnames <- colnames(data)
  data2 <- data[, c(subject, outcomes)]
  p <- length(outcomes)
  if( p < 2){
    data2$x_NA <- is.na(data2[,c(outcomes)])
    data2$x_NA <- as.numeric(data2$x_NA)
  }
  if(p>=2){
    data2$x_NA <- apply(data2[,c(-which(colnames %in% subject))],MARGIN = 1, FUN = AnyObsByLine)
  }
  data <- data[which(data2$x_NA ==0),]
  return(data)
}

#=====================================================================================
#' function that discretized a vector of time with a given delta.
#'
#' @param Time input continuous time vector
#' @param delta discretized step
#'
#' @return a discretized  time vector
#' 
f_TempsDiscr <- function(Time, Delta){
  AxeT <-seq(from =min(Time), to = max(Time), by = Delta)
  Time_d <- rep(0,length(Time))
  
  for(j in 1:length(Time)){
    i <-0
    Time_d[j] <- NA
    if(!is.na(Time[j])){
      Time_d[j] <- 0
      while(Time[j] > Delta*i){
        i <- i+1
      }
      Time_d[j] <- i*Delta
    }
  }
  return(Time_d)
}


#=====================================================================================
#' function that discretized a vector of time with a given delta.
#'
#' @param rdata input data which time point has to be discretized
#' @param outcomes a vector of outcomes names
#' @param predictors independent variables  to be be included in the modeling
#' @param subject subject identifiant
#' @param Time colname indiquating the time
#' @param Delta discretized time step 
#'
#' @return a discretized  time vector
#' 
DataDiscretization <-function(rdata, subject, outcomes, predictors = NULL, Time, Delta){
  
  cl <- match.call()
  colnames<-colnames(rdata)
  
  # Is subject colname available in the data
  if(!(subject %in% colnames))stop("Data discretisation failed: Subject colname should be in the data \n")
  
  # Is predictors colnames available in the data
  if(!(all(predictors %in% colnames)))stop("Data discretisation failed: All predictors colnames should be in the data \n")
  
  # Is Time colname available in the data
  if(!(unique(Time) %in% colnames))stop("Data discretisation failed: Time colname should be in the data \n")
  
  # # check if Delta is not greater than min(Time_(j+1)-Time(j))
  # Subjects = unique(rdata[,subject])
  # min_Time_diff = 1e23
  # for(n in 1:length(Subjects)){
  #   for(k in 1:length(outcomes)){
  #   Time_df = sort(rdata[which(rdata[,subject]==Subjects[n] & !is.na(rdata[,outcomes[k]])),Time])
  #   #
  #   for(i in 2:length(Time_df)){
  #     min_Time_diff <- min(min_Time_diff,  Time_df[i]-Time_df[(i-1)], na.rm=TRUE)
  #   }
  #   }
  # }
  # if(min_Time_diff > Delta) stop("Discretization failed: Discretization value could not be greater than the delay between two visit time")
  # 
  Time = rep(unique(Time),length(outcomes))## replicate Time
  ##  pre-processing of data: retaining lines with at least one observed outcome value
  data <- OneOutcomeAtLeast(rdata, subject= subject, outcomes = outcomes)
  cat(paste("Number of rows of the initial data is:", dim(rdata)[1]),"\n", paste("After removing lines without any observation of the outcomes of interest, the number of rows is:", dim(data)[1]),"\n")
  
  K <- length(outcomes)
  T_k <-NULL
  nameTks <- NULL
  for(k in 1:K) {
    nameTk <- paste(Time[k], k, sep = "_")
    nameTks <- c(nameTks,nameTk)
    #discretization of the time by subject: on each subject specific dataset
    T_k <- cbind(T_k,assign(nameTk, f_TempsDiscr(data[,Time[k]], Delta)))
  }
  T_k <- data.frame(T_k)
  colnames(T_k) <-nameTks
  data2 <-data.frame(cbind(data[,c(subject,outcomes, predictors)],T_k))
  
  # merge
  data3 <- na.omit(data2[,c(subject, outcomes[1], predictors, nameTks[1])])
  colname <- colnames(data3)
  # colname[which(colname==nameTks[1])] <- paste(Time[1],"d", sep = "_")
  colname[which(colname==nameTks[1])] <- Time[1]
  colnames(data3) <- colname
  
  if(K>1){
    for(k in 2:K){
      data_k <- na.omit(data2[,c(subject, outcomes[k], predictors, nameTks[k])])
      # changing the name of T_i colname into T for the merging
      colname <- colnames(data_k)
      # colname[which(colname==nameTks[k])] <- paste(Time[1],"d", sep = "_")
      colname[which(colname==nameTks[k])] <- Time[1]
      colnames(data_k) <- colname
      
      data3 <- merge(data3,data_k, by=c(subject,predictors,Time[1]), all.x = TRUE, all.y = TRUE)
      data3 <- data3[order(data3[,subject], data3[,Time[1]]), ]
    }
  }
  
  # At a visit, an outcome could be missing, but not all outcomes.
  # That is, for each row must have at least non missing observation.
  # So we need to check this after the discretization process.
  
  data4 <- OneOutcomeAtLeast(data=data3, subject=subject, outcomes = outcomes)
  if(dim(data4)[1] < dim(data3)[1])stop("Discretization failed: After discretization, All marker are missing at some visits")
  # check if Delta is not greater than min(Time_(j+1)-Time(j))
  if(dim(unique(data4))[1] != dim(data4)[1]) stop("Discretization failed: Some rows are the same in the dataset, because of a too large discretisation argument")
  
  return(data4)
}