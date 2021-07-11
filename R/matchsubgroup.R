
#'@title Re-screen the initial time slice of the subject population and data according to the conditions
#'@description The data of the subjects that meet the conditions for the first time and the data of
#'subsequent time slices will also be selected. According to the given screening conditions, the eligible subjects are screened
#'@usage matchsubgroup(
#'time,
#'status,
#'timeslicedata,
#'tspatientid,
#'varname,
#'varvalue
#')
#'@param time list object;Elements sorted by time node, each element is a
#'Dataframe object, representing event time or censoring time for subjects
#'@param status list object;Elements sorted by time node, each element is a
#'Dataframe object, representing status, 1 if death or event, 0 otherwise.
#'@param  timeslicedata list object; Elements sorted by time node, each element is a
#'Dataframe object, representing risk factors for the subject
#'@param tspatientid list object; Elements sorted by time node, each element represents the subjects identification number
#'@param  varname list object;The variable used to screen subjects, and the variable
#'needs to be a variable in timeslicedata
#'@param varvalue  list object;Subjects whose varname variable is equal to varvalue will be selected
#'@details According to the input time, status, variables, subject ID, etc., the data of eligible
#'subjects is screened through specified conditions, and the subject and the variable data of the
#'first and subsequent time slices are re-screened. The final returned result contains four list
#'objects: time, state, timeslicedata, subject ID(tspatientid).
#'@return Returns a list of the following items:
#'time,state,timeslicedata,tspatientid
#'@author Shen Lujun and ZhangTao
#'@export
#'@examples
#'varname=list('Amount of Hepatic Lesions','Largest Diameter of Hepatic Lesions (mm)') ,varvalue=list(1,1)
#'df <- matchsubgroup(time,status,tsdata,tspatientid,varname=
#'                       list('Amount of Hepatic Lesions') ,varvalue=list(1))
#'

matchsubgroup <- function(time,status,timeslicedata,tspatientid,varname,varvalue){

  d_time <- data.frame()
  d_status <- data.frame()
  d_tsdata <- data.frame()
  d_tsid <- data.frame()
  d_timenode <- data.frame()
  for (i in 1:length(time)) {
    timenode <- rep(i,dim(time[[i]])[1])
    timenode <- data.frame(timenode)
    d_timenode <- rbind(d_timenode,timenode)
    d_time <- rbind(d_time,time[[i]])
    d_status <- rbind(d_status,status[[i]])
    d_tsdata <- rbind(d_tsdata,timeslicedata[[i]])
    d_tsid <- rbind(d_tsid,tspatientid[[i]])
  }

  dataset <- cbind(d_time,d_status,d_tsdata,d_tsid,d_timenode)
  names(d_tsid) <- "id_id"
  dataset <- cbind(dataset,d_tsid)

  names(d_time) <- "time_time"
  dataset <- cbind(dataset, d_time)

  dataset <- dataset %>% group_by(id_id) %>% arrange(id_id,timenode)

  data <- dataset

  if (!is.list(varname)) {
    stop('varname should be a list!')
  }
  if (!is.list(varvalue)) {
    stop('varvalue should be a list!')
  }

  if (length(varname)!=length(varvalue)) {
    stop('The length ofvarvalue should equal with the length of varvalue!')
  }

  len <- length(varname)

  if (len==0){
    stop("Variables not specified for selecting data")
  }


  for (i in 1:len) {

    dataset <- dataset[which(dataset[,varname[[i]]]==varvalue[[i]]),]
  }

  ind <- duplicated(dataset$id_id)
  id_index <- subset(dataset,select = c("id_id","timenode"))
  id_index = data.frame(id_index)
  id_index = id_index[!ind,]

  n_data = data.frame()
  for (i in 1:dim(id_index)[1]) {
    m_data = data[which(data$id_id==id_index$id_id[i]),]
    m_data = m_data[which(m_data$timenode<=id_index$timenode[i]),]
    n_data = rbind(n_data,m_data)
  }

  timename = names(time[[1]])
  statusname = names(status[[1]])
  varname = names(timeslicedata[[1]])

  num <- length(time)

  time <- list()

  status <- list()

  timeslicedata <- list()

  tspatientid <- list()

  for (i in 1:num){

    data <- n_data[n_data['timenode']==i,]
    if(dim(data)[1]==0){
      break()
    }

    time <- c(time,list(data[timename]))

    status <- c(status,list(data[statusname]))

    tspatientid <- c(tspatientid,list(data['id_id']))

    c_data <- subset(data, select =varname)

    timeslicedata <- c(timeslicedata,list(c_data))
  }

  result <- list()
  result$time = time
  result$status = status
  result$tspatientid = tspatientid
  result$timeslicedata = timeslicedata
  return(result)
}
