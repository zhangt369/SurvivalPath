#  Reclassify variables into binary variables based on survivalROC
# Based on the survivalROC package, using the survival time
#and survival state, calculate the optimal cutoff for the binary classification
#of each variable using data at the first time slice. According to the cutoff, the variables of
#different nodes are classified into two categories.
# classifydata(
#time,
#status,
#timeslicedata,
#tspatientid,
#predict.time,
#isfill=TRUE
#)
#@param time Event time or censoring time for subjects
#@param status Indicator of status, 1 if death or event, 0 otherwise
#@param  timeslicedata list object; each dataframe of the list represent data at individual time slices
#@param tspatientid list object; The unique identification number corresponding to the subject at each time slice
#@param  predict.time Time point to predict during drawing the ROC curve
#@param isfill  Logical value, used to confirm whether to fill in missing data. If it is True, then fill
#@details According to the first time node, the survival data is used to calculate the optimal classification threshold
#of different risk factors. The risk factor threshold calculated at the first time node is used to calculate
#the grading of the risk factors at each time node. Return the cutoff corresponding to each risk factor and
#the classified timeslicedata.
#@return Returns a list of the following items:
#@return timeslicedata Risk factors after being classified
#@return cutoff  Optimal classification threshold
#@seealso survivalROC
#@import survivalROC



classifydata <- function(time,status,timeslicedata,tspatientid,predict.time,isfill=TRUE){

  varnames <- names(timeslicedata[[1]])

  variables <- timeslicedata[[1]]

  result <- list()
  cutoff <- list()
  cut.variable <- list()

  for (v in 1:length(varnames)){

    varn <- varnames[v]
    #print(varn)
    variable <- variables[varn]
    isnaindex <- is.na(variable)

    #print(names(variable))

    variable <- unlist(variable)

    num = length(unique(variable))

    if (num > 2){

      Stime = time[[1]][[1]][!isnaindex]
      Sstatus = status[[1]][[1]][!isnaindex]
      marker = variable[!isnaindex]

      #if(varn=="AFP"){
      #  thresholds=640.3
      #}else{

      roc= survivalROC(Stime=Stime,status=Sstatus, marker = marker,predict.time = predict.time,method="KM")
      ## roc <- roc(status[[1]][[1]], variable,main="Confidence intervals", percent=TRUE,ci=TRUE )
      ## ci <- ci(roc, of="thresholds", thresholds="best")
      ## thresholds <- attr(ci,"thresholds")[1]

      thresholds=roc$cut.values[which.max(abs(roc$TP-roc$FP))]
      #}
      #print(thresholds)

      for (i in 1:length(timeslicedata)){

        if (isfill & sum(is.na(timeslicedata[[i]][[v]]))>0) {

          timeslicedata[[i]][[v]][is.na(timeslicedata[[i]][[v]])] = mean(timeslicedata[[i]][[v]],na.rm=TRUE)

        }

        # data.frame


        timeslicedata[[i]][[v]] <- ifelse(timeslicedata[[i]][[v]] <= thresholds,0,1)
      }

      cutoff <- c(cutoff,thresholds)
      #print(cutoff)
      cut.variable <- c(cut.variable,varn)
    }else if(num==1){

      #print("single")

      thresholds = variable[1]
      for (i in 1:length(timeslicedata)){

        if (isfill & sum(is.na(timeslicedata[[i]][[v]]))>0) {

          timeslicedata[[i]][[v]][is.na(timeslicedata[[i]][[v]])] = mean(timeslicedata[[i]][[v]],na.rm=TRUE)

        }

        # data.frame

        timeslicedata[[i]][[v]] <- ifelse(timeslicedata[[i]][[v]] <= thresholds,0,1)
      }
      cutoff <- c(cutoff,thresholds)
      #print(cutoff)
      cut.variable <- c(cut.variable,varn)
    }

    cutoff <- data.frame(cutoff)
    names(cutoff)  <- cut.variable

  }

  (result$timeslicedata <- timeslicedata)
  (result$cutoff <- cutoff)
  result
}

