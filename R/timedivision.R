
#'@title  According to the specified time of the original data, select the data that conforms to the periodicity,
#'and generate the data format that conforms to the survivalpath function
#'@description Perform data preprocessing for survivalpath, filter and regularize
#'data entries that conform to continuous time intervals according to time.
#'@usage timedivision(
#'dataset,
#'ID,
#'time,
#'period=30,
#'left_interval = 0.5,
#'right_interval = 1.5
#')
#'@param dataset Refer to the built-in data set of this package
#'@param ID Subject ID corresponding to each row of data in the dataset
#'@param  time Represents the time when each row of data in the dataset is collected
#'@param period Represents the interval at which subjects with the same identification number are screened
#'@param left_interval The data collection time is less than the difference between the node time and
#'left_interval*period, that is, the value is considered to belong to the time node. For data collected
#'multiple times at the same time node, the earliest data is selected to represent the data collected by the node.
#'@param right_interval The data collection time is less than the difference between right_interval*period and the
#'node time, that is, the value is considered to belong to the time node. For data collected
#'multiple times at the same time node, the earliest data is selected to represent the data collected by the node.
#'@details for time point
#'@return dataframe;Data with the same ID number, sorted by consecutive time nodes
#'@author Shen Lujun and ZhangTao
#'@examples
#'data("dataset")
#'dataset = timedivision(X2021data,"ID","Date",period = 90,left_interval = 0.5,right_interval=1.5)
#'@export
#'
timedivision <- function(dataset,ID,time,period=30,left_interval = 0.5, right_interval = 1.5){
  if (left_interval>1. | left_interval<=0) {
    stop('left_interval should be between 0 and 1!')
  }
  if (right_interval<=1. | right_interval>=2) {
    stop('right_interval should be between 1 and 2!')
  }
  if ((left_interval+right_interval)<=1 |(left_interval+right_interval)>2) {
    stop('The sum of left_interval and right_interval should be between 1 and 2!')
  }

  dataset$time_time <- as.Date(dataset[[time]])
  dataset$id_id <- dataset[[ID]]

  # Access to multiple occurrences of data
  test <- dataset %>% group_by(id_id) %>% filter(n() > 1) %>% arrange(id_id,time_time) #6890

  # Gets data that appears only once
  test_1 <- dataset %>% group_by(id_id) %>% filter(n() == 1) %>% arrange(id_id,time_time) #85
  test_1$timenode = 1

  # According to the time and Id number, get the time and Id number
  #test1 <- dataset %>% group_by(id_id) %>% slice(1) %>% arrange(id_id,time_time) #1146

  #Remove Duplicates

  print("----------------------dup-----------------------")
  print(dim(test)[1])

  # time2: Number of cycles(period) between each time point and the minimum time point
  # time3: Difference between the number of intervals between each time point and the subsequent time point
  test1 <- test %>% mutate(time2 = as.numeric(time_time - min(time_time, na.rm = TRUE))/period)%>%
    mutate(time3 = ifelse(((time2 - lag(time2)) <= right_interval) & ((time2 - lag(time2)) >= left_interval),1,round(time2 - lag(time2))) ) %>%
    filter( is.na(time3) == TRUE | time3 ==1 ) %>%
    arrange(id_id,desc(time_time)) %>%
    ungroup

  # The last step is to re-group the selected data
  test <- test1 %>% group_by(id_id) %>% arrange(id_id,time_time)


  current_dim1 = 0
  print(dim(test)[1])
  print(current_dim1)
  while (current_dim1 != dim(test)[1]) {

    current_dim1 = dim(test)[1]

    test1 <- test %>% mutate(time4 = as.numeric(time_time - min(time_time, na.rm = TRUE))/period)%>%
      mutate(time5 = ifelse(((time4 - lag(time4)) <= right_interval) & ((time4 - lag(time4)) >= left_interval),1,round(time2 - lag(time2)) ) ) %>%
      filter( is.na(time5) == TRUE | time5 ==1 ) %>%
      arrange(id_id,desc(time_time)) %>%
      ungroup

    test <- test1 %>% group_by(id_id)  %>% arrange(id_id,time_time)
  }

  test1 <- test %>% group_by(id_id) %>% arrange(id_id,time_time) %>% mutate(time6 = 1) %>% mutate(timenode = cumsum(time6)) %>% ungroup

  test1 <- bind_rows(test_1,test1)

  test1 <- subset(test1,select = -c(time2,time3,time4,time5,time6,id_id,time_time))

  test1
}
