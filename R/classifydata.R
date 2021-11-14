#'@title  Reclassify variables into binary variables based on survivalROC
#'@description Based on the survivalROC package, using the survival time
#'and survival state, calculate the optimal cutoff for the binary classification
#'of each variable using data at the first time slice. According to the cutoff, the variables of
#' different nodes are classified into two categories.
#'@usage classifydata(
#'time,
#'status,
#'timeslicedata,
#'tspatientid,
#'predict.time,
#'isfill=T
#')
#'@param time Event time or censoring time for subjects
#'@param status Indicator of status, 1 if death or event, 0 otherwise
#'@param  timeslicedata list object; each dataframe of the list represent data at individual time slices
#'@param tspatientid list object; The unique identification number corresponding to the subject at each time slice
#'@param  predict.time Time point to predict during drawing the ROC curve
#'@param isfill  Logical value, used to confirm whether to fill in missing data. If it is True, then fill
#'@details According to the first time node, the survival data is used to calculate the optimal classification threshold
#'of different risk factors. The risk factor threshold calculated at the first time node is used to calculate
#'the grading of the risk factors at each time node. Return the cutoff corresponding to each risk factor and
#'the classified timeslicedata.
#'@return Returns a list of the following items:
#'@return timeslicedata Risk factors after being classified
#'@return cutoff  Optimal classification threshold
#'@seealso survivalROC
#'@export
#'@import survivalROC
#'@examples
#'data("DTSDHCC")
#'
#'dataset = timedivision(DTSDHCC,"ID","Date",period = 90,left_interval = 0.5,right_interval=0.5)
#'
#'time <- list()
#'
#'status <- list()
#'
#'tsdata <- list()
#'
#'tsid <- list()
#'
#'treatment <- list()
#'
#'for (i in 1:10){
#'
#'  data <- dataset[dataset['time_slice']==i,]
#'
#'  time <- c(time,list(data['OStime_day']))
#'
#'  status <- c(status,list(data['Status_of_death']))
#'
#'  tsid <- c(tsid,list(data['ID']))
#'
#'  c_data <- subset(data, select = c( "Age", "Amount of Hepatic Lesions", "Largest Diameter of Hepatic Lesions (mm)", "New Lesion",
#'    "Vascular Invasion" ,"Local Lymph Node Metastasis", "Distant Metastasis" , "Child_pugh_score" ,"AFP"))
#'
#'
#'  tsdata <- c(tsdata,list(c_data))
#'
#'  c_treatment <- subset(data, select = c("Resection"))
#'
#'  treatment <- c(treatment,list(c_treatment))
#'}
#'
#'
#' # predict.time
#'tsdata <- classifydata(time,status,tsdata,tsid,predict.time=365*1)
#'


classifydata <- function(time,status,timeslicedata,tspatientid,predict.time,isfill=T){

  varnames <- names(timeslicedata[[1]])

  variables <- timeslicedata[[1]]

  result <- list()
  cutoff <- list()
  cut.variable <- list()

  for (v in 1:length(varnames)){

    varn <- varnames[v]

    variable <- variables[varn]

    #print(names(variable))

    variable <- unlist(variable)

    num = length(unique(variable))

    if (num > 2){

      roc= survivalROC(Stime=time[[1]][[1]],status=status[[1]][[1]], marker = variable,predict.time = predict.time,method="KM")
      ## roc <- roc(status[[1]][[1]], variable,main="Confidence intervals", percent=TRUE,ci=TRUE )
      ## ci <- ci(roc, of="thresholds", thresholds="best")
      ## thresholds <- attr(ci,"thresholds")[1]

      thresholds=roc$cut.values[which.max(abs(roc$TP-roc$FP))]

      for (i in 1:length(timeslicedata)){

        if (isfill) {
          if (sum(is.na(timeslicedata[[i]][[v]]))>0){
            timeslicedata[[i]][[v]][is.na(timeslicedata[[i]][[v]])] = mean(timeslicedata[[i]][[v]],na.rm=T)
          }
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
