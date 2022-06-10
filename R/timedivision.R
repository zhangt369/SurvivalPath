
#'@title Convert Multiple Rows Arranged Time-Series Data into Time-Slices Data
#'@description Data preprocessing process essential for building survival path model. For each subject with  observations at different time point,
#'screen out specific observations at each specific time slice by setting associated parameters, includes period, left_interval and right_interval.
#'@usage timedivision(dataset,
#'ID,
#'time,
#'period=30,
#'left_interval = 0.5,
#'right_interval = 0.5
#')
#'@param dataset A multiple rows arranged time-series dataset, containing identification numbers, follow-up time points, risk factors, survival time, and survival status.
#'@param ID Character string, representing\code{ID} corresponding to each row of data in the dataset, which should be unique for each subject.
#'@param time Date format, which indicates time point of each observation.
#'@param period Numeric, utilized to customize follow-up sampling \code{period};normally counting in days.
#'@param left_interval Numeric, preferentially fall into the interval of (0,1). For a specific sampling in time slice \code{T}, the earliest sampling in the time interval [ \code{left_interval}*period, right_interval*period] is considered as the sampling data of the specific time slice \code{T}.
#'@param right_interval same as above.
#'@details This function is used to facilitate automatic generation of time-slice data. The date of observations for each subject should be arranged in ascending order.
#'The researchers can skip this process if they intend to prepare time-slice data manually or using customized codes. It's important to note that this
#'function only support data sampling of the "earliest" observation of interval in each time slice. If no observation fall into the interval of time
#'slice T, then sampling of observation in time slice T+1 for that subject will be terminated.
#'@return data.frame;observations of different time slices for each ID.The new data.frame returned added a new column "time_slice", which indicates the time slice of each observation included.
#'@export
#'@author Lujun Shen and Tao Zhang
#'@importFrom dplyr src
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


timedivision <- function(dataset,ID,time,period=30,left_interval = 0.5, right_interval = 0.5){

  dataset$time_time <- as.Date(dataset[[time]])
  dataset$id_id <- dataset[[ID]]
  ID <- unique(dataset$id_id)
  division <- function(x){
    data <- dataset[which(dataset$id_id==x),]
    data <- data%>% arrange(data$time_time)
    data$time2 <- as.numeric(data$time_time)
    num <- dim(data)[1]
    if (num==1){
      newdata <- data
      newdata$time_slice=1
      newdata
    }else if(num>1){
      index <- c(1)
      inittime <- data$time2[1]
      for (ii in 2:num) {
        left_interval_time <- inittime+period-left_interval*period
        right_interval_time <- inittime+right_interval*period+period
        currenttime <- data$time2[ii]
        if (currenttime>=left_interval_time & currenttime<= right_interval_time){
          index <- c(index,ii)
          inittime <- data$time2[1]+length(index)*period-period
          #print(paste(left_interval_time,currenttime,right_interval_time,sep = "-"))
        }
      }
      newdata <- data[index,]
      newdata$time_slice <- 1:length(index)
      newdata
    }
  }
  result <- lapply(ID, division)
  res <- bind_rows(result)
  #res <- subset(res,select=-c(id_id,time_time,time2))
  res

}
