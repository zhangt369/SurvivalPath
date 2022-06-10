#'@title Instantiate the an object of class Dynamic Time Series Data (DTSD)
#'@description Generate DTSD class objects using a dataframe. The dataframe should include unique identification number for each subject, multiple
#'rows arranged data (contain risk factors, survival time and outcomes) representing observations at different time slices/time points.
#'@usage generatorDTSD(dataset,
#'periodindex,
#'IDindex,
#'timeindex,
#'statusindex,
#'variable,
#'ifclassifydata=TRUE,
#'predict.time=365,
#'isfill=TRUE
#')
#'
#'@param dataset  A dataframe of time-series observations, containing identification numbers of each subject, index of time slice, value of risk factors, survival time, and survival outcomes.
#'@param periodindex Time slice indicator, represent index of time slice of specific observation, This variable is normally coded by integers, e.g. 0, 1, 2...
#'@param IDindex Variable name representing patient identification number.
#'@param timeindex  Variable name representing follow up time for censored data for each specific observation.
#'@param statusindex  The status indicator representing the patient's outcome status. For Overall survival, the status is normally coded by the policy 0=alive, 1=dead.
#'@param  variable List object containing the risk factors required for modeling.
#'@param  ifclassifydata A logical value, which is optional. Judgment on whether to classify risk factors automatically. When \code{ifclassifydata} is \code{TRUE} (default is TRUE), survivalROC method is used to find cutoff to dichotomize risk factors.
#'@param  predict.time Optional, Time of event assessment for identifying the best cutoff using survivalROC. When \code{ifclassifydata} is \code{TRUE}, \code{predict.time} is used in combination.
#'@param isfill  Logical value, used to confirm whether to fill in missing data. If it is True, then fill.
#'@return return a DTSD class object for survivalpath() function.
#'\item{time}{\code{time} list object; Event time or censoring time for subjects.Each element of the list represents, the event time or censoring time starting from each observation}
#'\item{status}{\code{status} list object; Indicator of status,normally use 0/1 coding. If death or event,1,otherwise,0. Each element of the list represents, the subject's outcome/event.}
#'\item{tsdata}{\code{tsdata} list object; Each element of \code{tsdata} contains the risk factors listed in \code{variable}. Each element of the list represents the data frame of each time slice, normally arranged in ascending order}
#'\item{tsid}{\code{tsid} list object; patient identification number. Each element of the list represents,the identification number of patient at each time slice}
#'\item{length}{  \code{time},\code{status},\code{tsdata} ,\code{tsid} are the same length \code{length}.}
#'\item{ts_size}{List object, representing sample size at each time slice.}
#'\item{cutoff}{ List object, representing cut-off values for each variable used for modeling.}
#'@details This function return a DTSD class object for conducting survivalpath function. This function facilitate enabling automatic binary classification of
#'continuous variables. When continuous variables need to be classified, survivalROC uses survival data at the predict.time to calculate cutoffs. The cutoff will
#'be used for construction of survival path at all time slices.
#'@examples
#'library(dplyr)
#'data("DTSDHCC")
#'id = DTSDHCC$ID[!duplicated(DTSDHCC$ID)]
#'set.seed(123)
#'id = sample(id,500)
#'miniDTSDHCC <- DTSDHCC[DTSDHCC$ID %in% id,]
#'dataset = timedivision(miniDTSDHCC,"ID","Date",period = 90,left_interval = 0.5,right_interval=0.5)
#'resu <- generatorDTSD(dataset,periodindex="time_slice",IDindex="ID" ,timeindex="OStime_day",
#'  statusindex="Status_of_death",variable =c( "Age", "Amount.of.Hepatic.Lesions",
#'  "Largest.Diameter.of.Hepatic.Lesions",
#'  "New.Lesion","Vascular.Invasion" ,"Local.Lymph.Node.Metastasis",
#'  "Distant.Metastasis" , "Child_pugh_score" ,"AFP"),predict.time=365*1)
#'
#'@export
#'



generatorDTSD <- function(dataset,periodindex,IDindex ,timeindex,
                 statusindex,variable ,
                 ifclassifydata=TRUE,predict.time=365,isfill=TRUE){

  time <- list()
  status <- list()
  tsdata <- list()
  tsid <- list()

  treatment <- list()
  timeslicesnum <- max(dataset$time_slice)
  for (i in 1:timeslicesnum){

    data <- dataset[dataset[periodindex]==i,]

    time <- c(time,list(data[timeindex]))

    status <- c(status,list(data[statusindex]))

    tsid <- c(tsid,list(data[IDindex]))

    c_data <- subset(data, select = variable)
    #c( "Age", "Amount of Hepatic Lesions", "Largest Diameter of Hepatic Lesions (mm)", "New Lesion",
    #"Vascular Invasion" ,"Local Lymph Node Metastasis", "Distant Metastasis" , "Child_pugh_score" ,"AFP")


    tsdata <- c(tsdata,list(c_data))

    c_treatment <- subset(data, select = c("Resection"))

    treatment <- c(treatment,list(c_treatment))
  }




  l_time <- length(time)
  l_status <- length(status)
  l_timesl <- length(tsdata)
  l_id <- length(tsid)


  ######################################################################################################
  #format

  if (!is.list(time)) {
    stop('time should be a list of data.frame!')
  }
  if (!is.list(status)) {
    stop('status should be a list of data.frame!')
  }
  if (!is.list(tsdata)) {
    stop('tsdata should be a list of data.frame!')
  }
  if (!is.list(tsid)) {
    stop('tsid should be a list of data.frame!')
  }

  #Data preprocessing, class 2 variables are transformed into class 2 variables

  #tsdata <- classifydata(time,status,tsdata,tsid)

  #Judging the number of input data time nodes is consistent
  if (l_time!=l_status | l_time!= l_timesl | l_time != l_id){
    stop(sprintf('The length of time,status,tsdata and tsid are:%d, %d, %d, %d,
                                but they should be the same length.',l_time,l_status,l_timesl,l_id))
  }



  #A patient's unique coding id required
  if (l_id==0){
    stop('The input can\'t be empty')
  }

  #A patient's unique coding id required,judgement
  patientID <- tsid[[1]]


  patientID_a <-unique(patientID[duplicated(patientID),])

  lenpid <- length(as.numeric(unlist(patientID_a)))
  if (lenpid!=0){
    stop(paste('The patient ID repeat. Every patient ID should be unique:',patientID_a))
  }



  len <- function(x){

    if("data.frame" %in%  class(x)){
      m= dim(x)[1]
    }
    if("list" %in% class(x)){
      m= length(x)
    }
    m
  }

  alllength <- unlist(lapply(tsid, len))

  if(ifclassifydata){

    # Reclassify variables into binary variables based on survivalROC
    tsdata_new <- classifydata(time,status,tsdata,tsid,predict.time=365*1,isfill=isfill)

    ndata <- list(time =time,
                  status = status ,
                  tsdata = tsdata_new$timeslicedata,
                  tsid = tsid,
                  length=length(tsid),
                  ts_size= alllength,
                  cutoff = tsdata_new$cutoff
    )

  }else{

    ndata <- list(time =time,
                  status = status ,
                  tsdata = tsdata,
                  tsid = tsid,
                  length=length(tsid),
                  ts_size= alllength,
                  cutoff = NULL
    )
  }


  #class(ndata) <- append(class(ndata),c("data","length","ts_size"))
  class(ndata) <- append(class(ndata),c("DTSD"))
  ndata

}


