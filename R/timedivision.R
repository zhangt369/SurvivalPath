

#'@title  Convert time-series data into time-slices in construction of survival paths
#'@description Perform data preprocessing for survivalpath, filter and regularize
#'data entries that conform to continuous time intervals (namely, into time slices) according to time.
#'@usage timedivision(
#'dataset,
#'ID,
#'time,
#'period=30,
#'left_interval = 0.5,
#'right_interval = 0.5
#')
#'@param dataset Refer to the colomn arranged time-series dataset, which has survival information
#'@param ID Subject ID corresponding to each row of data in the dataset, which should be unique for each participant
#'@param time the variable which indicates the time when each row of data in the dataset is collected
#'@param period The time interval at which time-series data of subjects with the same identification number are converted into time-slices data.
#'@param left_interval Refer to the times of time slice interval of lower level of time for inclusion the time point data into specific time slice.
#'For multiple time point data collected within the left and right interval of the specific time slice, the time point data closest to the integral
#'times of time slice interval and within the left interval is the first choice.
#'When there is no time point data in the left interval, the time point data closest to the integral times of time slice interval and within the right interval is the selected.
#'@param right_interval Refer to the times of time slice interval of higher level of time for inclusion the time point data into specific time slice.
#'When there is no time point data in the left interval, the time point data closest to the integral times of time slice interval and within the right interval is the selected.
#'@return dataframe;Data with the same ID number, sorted by consecutive time nodes
#'@author Lujun Shen and Tao Zhang
#'@import dplyr
#'@importFrom dplyr src
#'@examples
#'data("DTSDHCC")
#'dataset = timedivision(DTSDHCC,"ID","Date",period = 90,left_interval = 0.5,right_interval=0.5)
#'@export
#'
#'


timedivision <- function(dataset,ID,time,period=30,left_interval = 0.5, right_interval = 0.5){

  dataset$time_time <- as.Date(dataset[[time]])
  dataset$id_id <- dataset[[ID]]
  ID <- unique(dataset$id_id)
  division <- function(x){
    data <- dataset[which(dataset$id_id==x),]
    data <- data%>% arrange(time_time)
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

}
