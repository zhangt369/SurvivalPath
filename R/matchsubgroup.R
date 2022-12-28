
#'@title Screen and collect data of subjects that meet the given conditions
#'@description Based on screening criteria, for each specific subject, the data of the observation that meet the
#'conditions for the first time and the data of subsequent observations in following time slices will be collected.
#'The data from initial time slice that meet the given conditions to last time slice were then compiled into a new
#'time-slice dataset, with an aim to create personalized survival path map.
#'@usage matchsubgroup(
#'DTSD,
#'varname,
#'varvalue
#')
#'@param DTSD Object of class DTSD
#'@param  varname list object;The variable used to screen subjects, and the variables need to be contained in the time-slice data.
#'@param varvalue  list object;Subjects whose varname variable value equal to \code{varvalue} will be selected
#'@details According to the input time, status, variables, subject ID, etc., the data of eligible subjects is screened through specified
#'given conditions. The subject whose variable data of the first and subsequent time slices are sequentially screened. Once an observation
#'meet the given condition, data of that observation and the observations in following time slices will be for the subject will be collected.
#'Data of all subject that meet the criteria will be compiled into a new time-slice dataset. Based on the new dataset, the function returns a
#'new DTSD object was got. The final returned result contains four list objects: time, state, timeslicedata, subject ID (tspatientid).
#'@return Returns a DTSD object.
#'@author Shen Lujun and ZhangTao
#'@export
#'@examples
#'library(dplyr)
#'data("DTSDHCC")
#'id = DTSDHCC$ID[!duplicated(DTSDHCC$ID)]
#'set.seed(123)
#'id = sample(id,120)
#'miniDTSDHCC <- DTSDHCC[DTSDHCC$ID %in% id,]
#'dataset = timedivision(miniDTSDHCC,"ID","Date",period = 90,left_interval = 0.5,right_interval=0.5)
#'resu <- generatorDTSD(dataset,periodindex="time_slice",IDindex="ID" ,timeindex="OStime_day",
#'  statusindex="Status_of_death",variable =c( "Age", "Amount.of.Hepatic.Lesions",
#'  "Largest.Diameter.of.Hepatic.Lesions",
#'  "New.Lesion","Vascular.Invasion" ,"Local.Lymph.Node.Metastasis",
#'  "Distant.Metastasis" , "Child_pugh_score" ,"AFP"),predict.time=365*1)
#'
#'varname=list('Amount.of.Hepatic.Lesions')
#'varvalue=list(1)
#'df <- matchsubgroup(resu,varname=varname ,varvalue=varvalue)
#'
#'result <- survivalpath(df,time_slices =4)
#'

matchsubgroup <- function(DTSD,varname,varvalue){

  time <- DTSD$time
  status <- DTSD$status
  timeslicedata <- DTSD$tsdata
  tspatientid <- DTSD$tsid

  d_time <- data.frame()
  d_status <- data.frame()
  d_tsdata <- data.frame()
  d_tsid <- data.frame()
  d_timenode <- data.frame()
  for (i in 1:length(time)) {
    time_slice <- rep(i,dim(time[[i]])[1])
    time_slice <- data.frame(time_slice)
    d_timenode <- rbind(d_timenode,time_slice)
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

  dataset <- dataset %>% group_by(dataset$id_id) %>% arrange(dataset$id_id,time_slice)

  data <- dataset

  #if (!is.list(varname)) {
  #  stop('varname should be a list!')
  #}
  #if (!is.list(varvalue)) {
  #  stop('varvalue should be a list!')
  #}

  #if (length(varname)!=length(varvalue)) {
  #  stop('The length ofvarvalue should equal with the length of varvalue!')
  #}

  len <- length(varname)

  #if (len==0){
  #  stop("Variables not specified for selecting data")
  #}


  for (i in 1:len) {

    dataset <- dataset[which(dataset[,varname[[i]]]==varvalue[[i]]),]
  }

  ind <- duplicated(dataset$id_id)
  id_index <- subset(dataset,select = c("id_id","time_slice"))
  id_index = data.frame(id_index)
  id_index = id_index[!ind,]

  n_data = data.frame()
  for (i in 1:dim(id_index)[1]) {
    m_data = data[which(data$id_id==id_index$id_id[i]),]
    m_data = m_data[which(m_data$time_slice<=id_index$time_slice[i]),]
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

  alllength <- list()

  for (i in 1:num){

    data <- n_data[n_data['time_slice']==i,]
    if(dim(data)[1]==0){
      break()
    }

    time <- c(time,list(data[timename]))

    status <- c(status,list(data[statusname]))

    ID <- data['id_id']
    names(ID) = "ID"

    tspatientid <- c(tspatientid,list(ID))

    c_data <- subset(data, select =varname)

    timeslicedata <- c(timeslicedata,list(c_data))
    #print(dim(ID))
    alllength <- c(alllength,dim(ID)[1])
  }

  ndata <- list(time =time,
                status = status ,
                tsdata = timeslicedata,
                tsid = tspatientid,
                length=length(tspatientid),
                sublength= alllength,
                cutoff = DTSD$cutoff
  )
  class(ndata) <- append(class(ndata),c("DTSD"))

  return(ndata)
}
