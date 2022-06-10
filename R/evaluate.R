
#'@title Performance Evaluation of Survival Path Model
#'@description According to the survival path, using Harrell's concordance index C-index
#'to evaluate the discriminative ability of the survival path in each specified time slice for prognosis.
#'@usage evaluate(
#'survivalpath,
#'minnodesize
#')
#'
#'@param survivalpath  The output of the \code{survivalpath()}  function
#'@param minnodesize  The minimal sample size of specific node for inclusion of performance evaluation.
#'@return The evaluate function returns an object, which includes timeslice, Indexofmodes and Cindex
#'\item{timeslice}{time slice for which performance is evaluated.}
#'\item{Indexofnodes}{Nodes that is used for performance evaluation in the specified time slice.}
#'\item{Cindex}{Harrell's concordance index (C-index) value for specified time slices.}
#'@details For patients in each specific time slice, the class of survival path is regarded as a factor when computing the C-index value.
#'@importFrom Hmisc rcorr.cens
#'@importFrom stats predict
#'@export
#'@examples


evaluate <- function(survivalpath,minnodesize=10){
  dfresult = survivalpath$df

  maxpath <- survivalpath$maxpath

  #timeslice = timeslice+1

  res <- list()
  Cdx <- list()
  TS <- list()
  ION <- list()
  for (timeslice in 1:maxpath-1) {
    #print(timeslice)

    timeData = dfresult[dfresult$time_slice==timeslice & ! is.na(dfresult$sub_node),]
    tab <- table(timeData$sub_node)
    tab1 <- tab[tab>minnodesize]
    timeData <- timeData[timeData$sub_node %in%names(tab1),]
    if(dim(timeData)[1]==0){next()}
    if(dim(timeData)[2]==0){
      stop(paste("Time slice :",timeslice," does not exist in survival path",sep = ""))
    }

    if(length( unique(timeData$sub_node))==1){
      next()
    }

    timeData$node = as.factor(timeData$sub_node)

    treemodel <- rms::cph(Surv(timeData$OStime_day,timeData$Status_of_death)~`node`,data = timeData)

    risksocre <- predict(treemodel)

    cindex=1-rcorr.cens(risksocre,Surv(timeData$OStime_day,timeData$Status_of_death)) [[1]]

    Cdx <- c(Cdx,cindex)
    TS <- c(TS,timeslice)
    ION <- c(ION,list(unique(timeData$sub_node)))


  }

  res$timeslice <- TS
  res$Indexofnodes <- ION
  res$Cindex <- Cdx
  #print(res$Cindex)

  return(res)
}
